clear
% close all
rng shuffle
%% General settings

Pars.fc = 1e9; %carrier frequency
Ts=1/2*Pars.fc;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;%for URA generation

modOrder=8;
FFTLength=64;
NumSymbols=100;
snr = 3;

%initial position of UE 1 
x1o=-500;
%initial position of UE 2
x2o=0;
%step size for the movement of UE 1
step1=16.67;
%step size for the movement of UE 2
step2=0;

g.bs = [0;0;50]; %base station
g.t1 = [x1o;0;1.5]; %UE 1
g.t2 = [x2o ;0;1.5];%UE 2
x1=x1o;
x2=x2o;
for iter=1:40
    x1 = x1+step1; %do a step 
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
    Pars.fc = 1e9; %carrier frequency 
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
    chTaps=size(chan(1).delay); 
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
        chOutQUADRIGA=chOutQUADRIGA+chBS_T;%channel MP + interf
    end

    % no BF with QUADRIGA ch

    %add noise
    chOutQUADRIGA = awgn(chOutQUADRIGA, snr, 'measured');

    chOutQUADRIGA=transpose(chOutQUADRIGA);

    %store BER of UE1 no BF with quadriga ch
    bits=OFDMDemod(ofdmMod,chOutQUADRIGA(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1no_Q(iter)=(ratio);

    %store BER of UE2 no BF with quadriga ch
    bits=OFDMDemod(ofdmMod,chOutQUADRIGA(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber.UE2no_Q(iter)=(ratio);

    % LMS BF with QUADRIGA
    optimalWeight1 = LMSalgorithm(chOutQUADRIGA,waveform(:,1),numArrayElements);  
    optimalWeight2 = LMSalgorithm(chOutQUADRIGA,waveform(:,2),numArrayElements);   
    y1=chOutQUADRIGA*((optimalWeight1));
    y2=chOutQUADRIGA*((optimalWeight2));

    %store BER of UE1 LMS BF with quadriga ch
    bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1bf_Q(iter)=(ratio);

    %store BER of UE2 LMS BF with quadriga ch
    bits=OFDMDemod(ofdmMod,(y2),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber.UE2bf_Q(iter)=(ratio);

    %% free space channel

    %pathloss calculation for each terminal
    path_loss_t1 = ((4*pi*g.t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*g.t2_dist_BS)/Pars.lambda)^2;

    %compute rx signal at the BS and add noise
    receivedW = collectPlaneWave(g.BSarray, [waveform(:,1)*(1/sqrt(path_loss_t1)) waveform(:,2)*(1/sqrt(path_loss_t2))], [t1Angles' t2Angles'], Pars.fc);
    chOutFS = awgn(receivedW, snr, 'measured');

    %store BER of UE1 no BF with free space ch
    bits=OFDMDemod(ofdmMod,chOutFS(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1no_C(iter)=(ratio);

    %store BER of UE2 no BF with free space ch
    bits=OFDMDemod(ofdmMod,chOutFS(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber.UE2no_C(iter)=(ratio);

    %LMS BF
    optimalWeight1 = LMSalgorithm(chOutFS,waveform(:,1),numArrayElements);  
    optimalWeight2 = LMSalgorithm(chOutFS,waveform(:,2),numArrayElements);   
    y1=chOutFS*(optimalWeight1);
    y2=chOutFS*(optimalWeight2);


    %store BER of UE1 LMS BF with free space ch
    bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1bf_C(iter)=(ratio);

    %store BER of UE2 LMS BF with free space ch
    bits=OFDMDemod(ofdmMod,(y2),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    ber.UE2bf_C(iter)=(ratio);

end



ber.figureUE1_QvsC=figure('name','BER LMS quadriga vs free space UE1');
semilogy(linspace(x1o,x1,length(ber.UE1bf_Q)),ber.UE1bf_Q);

hold on
semilogy(linspace(x1o,x1,length(ber.UE1bf_C)),ber.UE1bf_C);

hold off
legend('quadriga','free sapce')
ylim([0.0000001 1])


ber.figureUE2_QvsC=figure('name','BER LMS quadriga vs free space UE2');
semilogy(linspace(x2o,x2,length(ber.UE2bf_Q)),ber.UE2bf_Q);
hold on
semilogy(linspace(x2o,x2,length(ber.UE2bf_C)),ber.UE2bf_C);
hold off
legend('quadriga','free space','location','best')
ylim([1e-5 1])
xlabel('Distance (BS in the origin)')
ylabel('BER')






