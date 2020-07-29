clear
close all
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
snr = 0;

%initial position  
x1o=-500;

%step size for the movement
step1=25;


g.bs = [0;0;50]; %base station position
x1=x1o;
g.t1 = [x1;0;1.5];

for iter=1:40
    x1 = x1+step1; %do a step 

    g.t1 = [x1;0;1.5];


    %azimuth angle
    g.az_t1=rad2deg(atan2(g.t1(1),g.t1(2)));


    %elevation angle
    g.el_t1=rad2deg(atan2(g.bs(3),sqrt(g.t1(1)^2+g.t1(2)^2)));


    %relevant angles for each terminal
    t1Angles = [g.az_t1 g.el_t1];


    %distances from base station in xy plane (BS is in (0,0))
    g.t1_dist = sqrt(g.t1(1)^2+g.t1(2)^2);

    %signal path
    g.t1_dist_BS=sqrt(g.bs(3)^2+g.t1_dist^2);


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
    l.no_tx = 1;
    l.rx_track.name= 'BS';
    l.rx_position = g.bs; %position fetching

    tx_track1= qd_track('linear',0,0); %track of the 1st terminal with random length between 0 and 20 m and random direction in [0,2*pi]

    tx_track1.name = 'UE1';

    l.tx_track(1,1) = copy(tx_track1); %track establishment


    l.tx_position=[g.t1];
    %  pos=l.visualize();
    pos.Name='Positions';
    pos.NumberTitle='off';
    l.set_pairing; %pair all tx/rx links
    chan = l.get_channels(); %channel generation


    %OFDM tx signals
    [ofdmMod,waveform,in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder);


    %OFDM demodulators
    title='TX Constellation by terminals';
    OFDMDemod(ofdmMod,waveform,modOrder,false,title);

    % channels convolution
    chTaps=size(chan(1).delay); 
    chBS_T=zeros(chTaps(1),length(waveform));
    TsVect=0:Ts:Ts*(length(waveform)-1);

    for antenna=1:1:chTaps(1)
        
        for path=1:1:chTaps(3)
            
            inX=TsVect-chan(1).delay(antenna,1,path,1);
            inY=interp1(TsVect,waveform,inX,'pship');
            chBS_T(antenna,:)=inY*chan(1).coeff(antenna,1,path,1)+chBS_T(antenna,:);%channel bs-terminal

        end
        
    end

    % no BF QUADRIGA

    %add noise
    chOutQUADRIGA = awgn(chBS_T, snr, 'measured');

    chOutQUADRIGA=transpose(chOutQUADRIGA);
    
    %store BER no BF with quadriga ch
    bits=OFDMDemod(ofdmMod,chOutQUADRIGA(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1no_Q(iter)=(ratio);

    %  LMS algorithm with QUADRIGA ch
    optimalWeight = LMSalgorithm(chOutQUADRIGA,waveform,numArrayElements);  
    y1=chOutQUADRIGA*((optimalWeight));
    
    %store BER with LMS BF and quadriga ch
    bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1bf_Q(iter)=(ratio);



    %% free space channel

    %pathloss calculation 
    path_loss_t1 = ((4*pi*g.t1_dist_BS)/Pars.lambda)^2;
    
    %compute rx signal at the BS and add noise
     receivedW = collectPlaneWave(g.BSarray, [waveform*(1/sqrt(path_loss_t1))], [t1Angles'], Pars.fc);
     chOutFS = awgn(receivedW, snr, 'measured');

    %store BER no BF with free space ch
    bits=OFDMDemod(ofdmMod,chOutFS(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1no_C(iter)=(ratio);



    % LMS BF
    optimalWeight = LMSalgorithm(chOutFS,waveform,numArrayElements);  
    y1=chOutFS*(optimalWeight);



    %store BER LMS BF with free space ch
    bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    ber.UE1bf_C(iter)=(ratio);


end

% plot the BER as a function of the distance in semilog scale
ber.figureUE1_QvsC=figure('name','BER LMS quadriga vs free space single UE');
semilogy(linspace(x1o,x1,length(ber.UE1bf_Q)),ber.UE1bf_Q);
hold on
semilogy(linspace(x1o,x1,length(ber.UE1bf_C)),ber.UE1bf_C);
hold off
legend('quadriga','free space')
ylim([1e-10 1])
xlabel('Distance (BS in the origin)')
ylabel('BER')







