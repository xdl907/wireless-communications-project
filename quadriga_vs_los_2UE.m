clear
% close all
rng shuffle
%% General settings

Pars.fc = 1e9; %carrier frequency
Ts=1/2*Pars.fc;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;

modOrder=8;
FFTLength=64;
NumSymbols=100;
snr = 3;

x1o=-700;
x2o=0;
step1=17.5;
step2=0;

g.bs = [0;0;50]; %base station
g.t1 = [x1o;0;1.5]; %terminal 1
g.t2 = [x2o ;0;1.5];
x1=x1o;
x2=x2o;
for iter=1:40
% g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal 1
% g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal 2
x1 = x1+step1; %terminal 1
x2 = x2+step2;
g.t1 = [x1;0;1.5];
g.t2 = [x2;0;1.5];

%azimuth angles
g.az_t1=rad2deg(atan2(g.t1(1),g.t1(2)));
g.az_t2=rad2deg(atan2(g.t2(1),g.t2(2)));


%elevation angles
g.el_t1=rad2deg(atan2(g.bs(3),sqrt(g.t1(1)^2+g.t1(2)^2)));
g.el_t2=rad2deg(atan2(g.bs(3),sqrt(g.t2(1)^2+g.t2(2)^2)));


%relevant angles for each terminal
t1Angles = [g.az_t1 g.el_t1];
t2Angles = [g.az_t2 g.el_t2];


%distances from base station in xy plane (BS is in (0,0))
g.t1_dist = sqrt(g.t1(1)^2+g.t1(2)^2);
g.t2_dist = sqrt(g.t2(1)^2+g.t2(2)^2);


%signal path
g.t1_dist_BS=sqrt(g.bs(3)^2+g.t1_dist^2);
g.t2_dist_BS=sqrt(g.bs(3)^2+g.t2_dist^2);


%% Phased array definition
Pars.fc = 1e9; %carrier frequency (add LTE one?)
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc; %wavelength 
g.BSarray = phased.URA('Size',[4 4], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x'); %URA definition
g.bs_antenna_pos=getElementPosition(g.BSarray); %positions of the antennas in the array

%% Quadriga channel model
l = qd_layout;
l.set_scenario('QuaDRiGa_UD2D_LOS'); %propagation model

array_tx = qd_arrayant('omni');
array_rx = qd_arrayant('omni');
array_rx.no_elements = 16;
array_rx.element_position=g.bs_antenna_pos; %fetching the URA antenna positions for the quadriga layout

l.tx_array = array_tx; %array definitions
l.rx_array = array_rx;
l.no_rx = 1;
l.no_tx = 2;
l.rx_track.name= 'BS';
l.rx_position = g.bs; %position fetching

tx_track1= qd_track('linear',0,0); %track of the 1st terminal with random length between 0 and 20 m and random direction in [0,2*pi]
tx_track2= qd_track('linear',0,0); %track of the 2nd terminal with random length between 0 and 30 m and random direction in [0,2*pi]


tx_track1.name = 'UE1';
tx_track2.name = 'UE2';


l.tx_track(1,1) = copy(tx_track1); %track establishment
l.tx_track(1,2) = copy(tx_track2);

l.tx_position=[g.t1,g.t2];
%  pos=l.visualize();
pos.Name='Positions';
pos.NumberTitle='off';
l.set_pairing; %pair all tx/rx links
chan = l.get_channels(); %channel generation


%OFDM tx signals

[ofdmMod,waveform(:,1),in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder);
[~,waveform(:,2),in_t2]=OFDMsignal(FFTLength, NumSymbols,modOrder);

%OFDM demodulators
title='TX Constellation by terminals';
OFDMDemod(ofdmMod,waveform(:,1),modOrder,false,title);

% channels convolution
chTaps=size(chan(1).delay); % channel 1&2 -> interf_Bs    channel 3 -> UE1_Bs    channel 4 -> UE2_Bs   
chOutQUADRIGA=zeros(chTaps(1),length(waveform(:,1)));
chBS_T0=zeros(chTaps(1),length(waveform(:,1)));
TsVect=0:Ts:Ts*(length(waveform(:,1))-1);

for channel=1:1:2
    chBS_T=chBS_T0; %reset
    for antenna=1:1:chTaps(1)
        for path=1:1:chTaps(3)
            
            inX=TsVect-chan(channel).delay(antenna,1,path,1);
            inY=interp1(TsVect,waveform(:,channel),inX,'pship');
            chBS_T(antenna,:)=inY*chan(channel).coeff(antenna,1,path,1)+chBS_T(antenna,:);%channel bs-terminal (only MP, no interf)
 
        end
    end
    chOutQUADRIGA=chOutQUADRIGA+chBS_T;%channel MP+interf (sum all obtained channels in one)
end

% no BF QUADRIGA

%add noise
chOutQUADRIGA = awgn(chOutQUADRIGA, snr, 'measured');

chOutQUADRIGA=transpose(chOutQUADRIGA);

title='Constellation without beamforming and QUADRIGA ch on UE1';
bits=OFDMDemod(ofdmMod,chOutQUADRIGA(:,end),modOrder,false,title);
[numbError,ratio]=biterr(in_t1,bits);
% fprintf('\nNumber of errors no beam 1 QUADRIGA ch: %d',numbError);
% fprintf('\nBER no beam 1 QUADRIGA ch: %f',ratio);
ber.UE1no_Q(iter)=(ratio);
% 
title='Constellation without beamforming and QUADRIGA ch on UE2';
bits=OFDMDemod(ofdmMod,chOutQUADRIGA(:,end),modOrder,false,title);
[numbError,ratio]=biterr(in_t2,bits);
% fprintf('\nNumber of errors no beam 2 QUADRIGA ch: %d',numbError);
% fprintf('\nBER no beam 2 QUADRIGA ch: %f',ratio);
ber.UE2no_Q(iter)=(ratio);

% algoritmo LMS QUADRIGA

optimalWeight1 = LMSalgorithm(chOutQUADRIGA,waveform(:,1),numArrayElements,Pars.lambda);  
 optimalWeight2 = LMSalgorithm(chOutQUADRIGA,waveform(:,2),numArrayElements,Pars.lambda);   

for i=1:length(waveform(:,1))
    y1(i)=chOutQUADRIGA(i,:)*((optimalWeight1));
     y2(i)=chOutQUADRIGA(i,:)*((optimalWeight2));
end

title='Constellation with LMS and QUADRIGA ch on UE1';
bits=OFDMDemod(ofdmMod,transpose(y1),modOrder,false,title);
[numbError,ratio]=biterr(in_t1,bits);
% fprintf('\nNumber of errors LMS1 QUADRIGA ch: %d',numbError);
% fprintf('\nBER LMS1 QUADRIGA ch: %f',ratio);
ber.UE1bf_Q(iter)=(ratio);

title='Constellation with LMS and QUADRIGA ch on UE2';
bits=OFDMDemod(ofdmMod,transpose(y2),modOrder,false,title);
[numbError,ratio]=biterr(in_t2,bits);
% fprintf('\nNumber of errors LMS2 QUADRIGA ch: %d',numbError);
% fprintf('\nBER LMS2 QUADRIGA ch: %f\n',ratio);
ber.UE2bf_Q(iter)=(ratio);

%% LOS channel

%pathloss calculation for each terminal
path_loss_t1 = ((4*pi*g.t1_dist_BS)/Pars.lambda)^2;
path_loss_t2 = ((4*pi*g.t2_dist_BS)/Pars.lambda)^2;

% 
% path_loss_t1 = 1;
% path_loss_t2 = 1;
% path_loss_i1 = 1;
% path_loss_i2 = 1;

  receivedW = collectPlaneWave(g.BSarray, [waveform(:,1)*(1/sqrt(path_loss_t1)) waveform(:,2)*(1/sqrt(path_loss_t2))], [t1Angles' t2Angles'], Pars.fc);
 chOutLOS = awgn(receivedW, snr, 'measured');
 
 title='Constellation without beamforming and LOS ch on UE1';
bits=OFDMDemod(ofdmMod,chOutLOS(:,end),modOrder,false,title);
[numbError,ratio]=biterr(in_t1,bits);
% fprintf('\nNumber of errors no beam 1 LOS ch: %d',numbError);
% fprintf('\nBER no beam 1 LOS ch: %f',ratio);
ber.UE1no_C(iter)=(ratio);


title='Constellation without beamforming and LOS ch on UE2';
bits=OFDMDemod(ofdmMod,chOutLOS(:,end),modOrder,false,title);
[numbError,ratio]=biterr(in_t2,bits);
%fprintf('\nNumber of errors no beam 2 LOS ch: %d',numbError);
%fprintf('\nBER no beam 2 LOS ch: %f',ratio);
ber.UE2no_C(iter)=(ratio);

% algoritmo LMS LOS

optimalWeight1 = LMSalgorithm(chOutLOS,waveform(:,1),numArrayElements,Pars.lambda);  
 optimalWeight2 = LMSalgorithm(chOutLOS,waveform(:,2),numArrayElements,Pars.lambda);   

for i=1:length(waveform(:,1))
    y1(i)=chOutLOS(i,:)*(optimalWeight1);
    y2(i)=chOutLOS(i,:)*(optimalWeight2);
end

title='Constellation with LMS and LOS ch on UE1';
bits=OFDMDemod(ofdmMod,transpose(y1),modOrder,false,title);
[numbError,ratio]=biterr(in_t1,bits);
% fprintf('\nNumber of errors LMS1 LOS ch: %d',numbError);
% fprintf('\nBER LMS1 LOS ch: %f',ratio);
ber.UE1bf_C(iter)=(ratio);
% 
title='Constellation with LMS and LOS ch on UE2';
bits=OFDMDemod(ofdmMod,transpose(y2),modOrder,false,title);
[numbError,ratio]=biterr(in_t2,bits);
%fprintf('\nNumber of errors LMS2 LOS ch: %d',numbError);
%fprintf('\nBER LMS2 LOS ch: %f\n',ratio);
ber.UE2bf_C(iter)=(ratio);

end


% ber.figureUE1no_C=figure('name','BER no BF LOS UE1');
% semilogy([1:1:iter],ber.UE1no_C);
% ylim([0.001 1])
% 
% ber.figureUE2no_C=figure('name','BER no BF LOS UE2');
% semilogy([1:1:iter],ber.UE2no_C);
% ylim([0.001 1])

% ber.figureUE1bf_C=figure('name','BER LMS LOS UE1');
% semilogy([1:1:iter],ber.UE1bf_C);
% ylim([0.000001 1])
% 
% ber.figureUE2bf_C=figure('name','BER LMS LOS UE2');
% semilogy([1:1:iter],ber.UE2bf_C);
% ylim([0.000001 1])


% ber.figureUE1no_Q=figure('name','BER no BF quadriga UE1');
% semilogy([1:1:iter],ber.UE1no_Q);
% ylim([0.0000001 1])

% ber.figureUE2no_Q=figure('name','BER no BF quadriga UE2');
% semilogy([1:1:iter],ber.UE2no_Q);
% ylim([0.0000001 1])
% 
% ber.figureUE1bf_Q=figure('name','BER LMS quadriga UE1');
% semilogy([1:1:iter],ber.UE1bf_Q);
% ylim([0.000001 1])
% 
% ber.figureUE2bf_Q=figure('name','BER LMS quadriga UE2');
% semilogy([1:1:iter],ber.UE2bf_Q);
% ylim([0.000001 1])


ber.figureUE1_QvsC=figure('name','BER LMS quadriga vs LOS UE1');
semilogy(linspace(x1o,x1,length(ber.UE1bf_Q)),ber.UE1bf_Q);
% set ( gca, 'xdir', 'reverse' )
hold on
semilogy(linspace(x1o,x1,length(ber.UE1bf_C)),ber.UE1bf_C);
% set ( gca, 'xdir', 'reverse' )
hold off
legend('quadriga','los')
ylim([0.0000001 1])
% set ( gca, 'xdir', 'reverse' )

ber.figureUE2_QvsC=figure('name','BER LMS quadriga vs LOS UE2');
semilogy(linspace(x2o,x2,length(ber.UE2bf_Q)),ber.UE2bf_Q);
hold on
semilogy(linspace(x2o,x2,length(ber.UE2bf_C)),ber.UE2bf_C);
hold off
legend('quadriga','los')
ylim([0.0000001 1])






