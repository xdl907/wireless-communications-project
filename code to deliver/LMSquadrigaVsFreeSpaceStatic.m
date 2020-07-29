clear;
%% general settings 
% Random coordinate generation
g.bs = [0,0,50]; %base station
g.t1 = [-50 + rand(2,1)*300;1.5]; %terminal 1
g.t2 = [-50 + rand(2,1)*300;1.5]; %terminal 2
g.i1 = [-50 + rand(2,1)*300;1.5]; %interferers 
g.i2 = [-50 + rand(2,1)*300;1.5];


Pars.fc = 1e9; %carrier frequency
Ts=1/2*Pars.fc;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;
%ofdm parameters
modOrder=8;
FFTLength=64;
NumSymbols=100; %can be increased for better precision in ber calculation

graph = false; %set to true to plot constellations, not really needed
%vector needed to store values of BER
ber1_nobf = [];
ber1_lms = [];
ber2_nobf = [];
ber2_lms = [];
ber3_nobf = [];
ber3_lms = [];
ber4_nobf = [];
ber4_lms = [];
snr_vec = [];

%OFDM tx signals
[~,waveform(:,1),in_i]=OFDMsignal(FFTLength, NumSymbols,modOrder);%signal tx by interf. 1
waveform(:,2)=waveform(:,1);%signal tx by interf. 2
[ofdmMod,waveform(:,3),in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder);%signal tx by UE 1
[~,waveform(:,4),in_t2]=OFDMsignal(FFTLength, NumSymbols,modOrder);%signal tx by UE 2

% waveform=[waveform_i,waveform_i,waveform_t1,waveform_t2];

% Spectrum Analyzers
% spectrum_t1 = dsp.SpectrumAnalyzer('SampleRate', 2*Pars.fc);
% spectrum_t1(waveform_t1);
% release(spectrum_t1);
% 
% spectrum_t2 = dsp.SpectrumAnalyzer('SampleRate', 2*Pars.fc);
% spectrum_t2(waveform_t2);
% release(spectrum_t2);
% 
% spectrum_i = dsp.SpectrumAnalyzer('SampleRate', 2*Pars.fc);
% spectrum_i(waveform_i);
% release(spectrum_i);
% 
% % OFDM Subcarrier Mapping
% showResourceMapping(ofdmMod);

%OFDM demodulators
title='TX Constellation by terminals';
bits=OFDMDemod(ofdmMod,waveform(:,3),modOrder,graph,title);


% MIMO array
g.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');
g.bs_antenna_pos=getElementPosition(g.BSarray); %positions of the antennas in the array

%% Array definitions for the Quadriga channel model
l = qd_layout;
l.set_scenario('QuaDRiGa_UD2D_LOS'); %propagation model

array_tx = qd_arrayant('omni');
array_rx = qd_arrayant('omni');
array_rx.no_elements = numArrayElements*numArrayElements;
array_rx.element_position=g.bs_antenna_pos; %fetching the URA antenna positions for the quadriga layout
l.tx_array = array_tx; %array definitions
l.rx_array = array_rx;
l.no_rx = 1;%n° of rx (base station)
l.no_tx = 4;%n° of tx
l.rx_track.name= 'BS';
l.rx_position = g.bs'; %position fetching
%% Quadriga channel model

tx_track1= qd_track('linear',20*rand,rand*2*pi); %track of the 1st terminal with random length between 0 and 20 m and random direction in [0,2*pi]
tx_track2= qd_track('linear',30*rand,rand*2*pi); %track of the 2nd terminal with random length between 0 and 30 m and random direction in [0,2*pi]
tx_track3= qd_track('linear',10*rand,rand*2*pi); %track of the 1st interf with random length between 0 and 10 m and random direction in [0,2*pi]
tx_track4= qd_track('linear',15*rand,rand*2*pi); %track of the 2nd interf with random length between 0 and 15 m and random direction in [0,2*pi]

tx_track1.name = 'UE1';
tx_track2.name = 'UE2';
tx_track3.name = 'INTERF 1';
tx_track4.name = 'INTERF 2';

l.tx_track(1,1) = copy(tx_track1); %track establishment
l.tx_track(1,2) = copy(tx_track2);
l.tx_track(1,3) = copy(tx_track3); 
l.tx_track(1,4) = copy(tx_track4);

l.tx_position=[g.t1,g.t2,g.i1,g.i2];
pos=l.visualize();
pos.Name='Positions';
pos.NumberTitle='off';
l.set_pairing; %pair all tx/rx links
chan = l.get_channels(); %channel generation


%azimuth angles
g.az_t1=rad2deg(atan2(g.t1(1),g.t1(2)));
g.az_t2=rad2deg(atan2(g.t2(1),g.t2(2)));
g.az_i1=rad2deg(atan2(g.i1(1),g.i1(2)));
g.az_i2=rad2deg(atan2(g.i2(1),g.i2(2)));

%elevation angles
g.el_t1=rad2deg(atan2(g.bs(3),sqrt(g.t1(1)^2+g.t1(2)^2)));
g.el_t2=rad2deg(atan2(g.bs(3),sqrt(g.t2(1)^2+g.t2(2)^2)));
g.el_i1=rad2deg(atan2(g.bs(3),sqrt(g.i1(1)^2+g.i1(2)^2)));
g.el_i2=rad2deg(atan2(g.bs(3),sqrt(g.i2(1)^2+g.i2(2)^2)));


%relevant angles for each terminal
t1Angles = [g.az_t1 g.el_t1];
t2Angles = [g.az_t2 g.el_t2];
i1Angles = [g.az_i1 g.el_i1];
i2Angles = [g.az_i2 g.el_i2];

%distances from base station in xy plane (BS is in (0,0))
g.t1_dist = sqrt(g.t1(1)^2+g.t1(2)^2);
g.t2_dist = sqrt(g.t2(1)^2+g.t2(2)^2);
g.i1_dist = sqrt(g.i1(1)^2+g.i1(2)^2);
g.i2_dist = sqrt(g.i2(1)^2+g.i2(2)^2);

%distances travelled by the signal
g.t1_dist_BS=sqrt(g.bs(3)^2+g.t1_dist^2);
g.t2_dist_BS=sqrt(g.bs(3)^2+g.t2_dist^2);
g.i1_dist_BS=sqrt(g.bs(3)^2+g.i1_dist^2);
g.i2_dist_BS=sqrt(g.bs(3)^2+g.i2_dist^2);


% channels convolution
for snr=5:15 %calculates channel with noise added for each SNR value in the interval 5-15
    chTaps=size(chan(3).delay); % channel 1&2 -> interf_Bs    channel 3 -> UE1_Bs    channel 4 -> UE2_Bs   
    chOutQ=zeros(chTaps(1),length(waveform(:,1)));
    chBS_T0=zeros(chTaps(1),length(waveform(:,1)));
    TsVect=0:Ts:Ts*(length(waveform(:,1))-1);

    for channel=1:1:4
        chBS_T=chBS_T0; %reset
        for antenna=1:1:chTaps(1)
            for path=1:1:chTaps(3)

                inX=TsVect-chan(channel).delay(antenna,1,path,1);
                inY=interp1(TsVect,waveform(:,channel),inX,'pship'); %interpolation with the tx signal
                chBS_T(antenna,:)=inY*chan(channel).coeff(antenna,1,path,1)+chBS_T(antenna,:);%coefficient multiplication: we get channel bs-terminal (only MP, no interference)

            end
        end
        chOutQ=chOutQ+chBS_T;%channel MP+interf (i.e. sum all obtained channels in one)
    end


    %add noise
    chOutQ = awgn(chOutQ, snr, 'measured');

    chOutQ=transpose(chOutQ);
    
    % no BF
    
    %Compute BER without BF on UE1
    title='Constellation without beamforming on UE1';
    bits=OFDMDemod(ofdmMod,chOutQ(:,end),modOrder,graph,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber1_nobf = [ber1_nobf, ratio]

    %Compute BER without BF on UE2
    title='Constellation without beamforming on UE2';
    bits=OFDMDemod(ofdmMod,chOutQ(:,end),modOrder,graph,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber2_nobf = [ber2_nobf, ratio]

    % LMS beamforming with quadriga ch

     %compute LMS weights
    optimalWeight1 = LMSalgorithm(chOutQ,waveform(:,3),numArrayElements);  
    optimalWeight2 = LMSalgorithm(chOutQ,waveform(:,4),numArrayElements); 
    
    %multiply rx signal by weights 
    y1=chOutQ*((optimalWeight1));
    y2=chOutQ*((optimalWeight2));

    title='Constellation with LMS on UE1';
    bits=OFDMDemod(ofdmMod,(y1),modOrder,graph,title);
    [numbError,ratio1]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
    fprintf('\nBER LMS1: %f',ratio1);
    ber1_lms = [ber1_lms, ratio1]%we append the new calculated values to the BER vector to plot them in the graphs

    title='Constellation with LMS on UE2';
    bits=OFDMDemod(ofdmMod,(y2),modOrder,graph,title);
    [numbError,ratio2]=biterr(in_t2,bits);
    fprintf('\nNumber of errors LMS2: %d',numbError);
    fprintf('\nBER LMS2: %f\n',ratio2);
    ber2_lms = [ber2_lms, ratio2]%we append the new calculated values to the BER vector to plot them in the graphs

    snr_vec = [snr_vec, snr]
    
    
    
    %% LOS channel

    %pathloss for each terminal
    path_loss_t1 = ((4*pi*g.t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*g.t2_dist_BS)/Pars.lambda)^2;
    path_loss_i1 = ((4*pi*g.i1_dist_BS)/Pars.lambda)^2;
    path_loss_i2 = ((4*pi*g.i2_dist_BS)/Pars.lambda)^2;


    receivedW = collectPlaneWave(g.BSarray, [waveform(:,3)*(1/sqrt(path_loss_t1)) waveform(:,4)*(1/sqrt(path_loss_t2)) waveform(:,1)*(1/sqrt(path_loss_i1)) waveform(:,2)*(1/sqrt(path_loss_i2))], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);
    chOutLOS = awgn(receivedW, snr, 'measured');
    
    % no BF
    
        %Compute BER without BF on UE1
    title='Constellation without beamforming and LOS ch on UE1';
    bits=OFDMDemod(ofdmMod,chOutLOS(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber3_nobf = [ber3_nobf, ratio]

    %Compute BER without BF on UE2
    title='Constellation without beamforming and LOS ch on UE2';
    bits=OFDMDemod(ofdmMod,chOutLOS(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber4_nobf = [ber4_nobf, ratio]

    % LMS BF with free space ch
    optimalWeight1 = LMSalgorithm(chOutLOS,waveform(:,3),numArrayElements);  
    optimalWeight2 = LMSalgorithm(chOutLOS,waveform(:,4),numArrayElements);   
   
    %multiply rx signal by weights
     y1=chOutLOS*((optimalWeight1));
     y2=chOutLOS*((optimalWeight2));

    title='Constellation with LMS and LOS ch on UE1';
    bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
     ber3_lms = [ber3_lms, ratio]
     
    title='Constellation with LMS and LOS ch on UE2';
    bits=OFDMDemod(ofdmMod,(y2),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber4_lms = [ber4_lms, ratio]

    
    
end

figure
ylim([0 1])
%plot the BER as a function of SNR in semilog scale
semilogy(snr_vec,ber1_lms,'-*',snr_vec,ber1_nobf,'-o',snr_vec,ber2_lms,'-*',snr_vec,ber2_nobf,'-o',snr_vec,ber3_lms,':*',snr_vec,ber3_nobf,':o',snr_vec,ber4_lms,':*',snr_vec,ber4_nobf,':o')
xlabel('SNR (dB)')
ylabel('BER')
legend('BER of UE1 with LMS BF and Quadriga','BER of UE1 without LMS BF and Quadriga','BER of UE2 with LMS BF and Quadriga','BER of UE2 without LMS BF and Quadriga',     'BER of UE1 with LMS BF and Free Space','BER of UE1 without LMS BF and Free Space','BER of UE2 with LMS BF and Free Space','BER of UE2 without LMS BF and Free Space')








