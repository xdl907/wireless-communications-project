clearvars
close all

%% general settings 
%random coordinate generation
bs = [0,0,50]; %base station
t1 = [-50 + rand(2,1)*100;0]; %terminal 1
t2 = [-50 + rand(2,1)*100;0]; %terminal 2
i1 = [-50 + rand(2,1)*100;0]; %interferers
i2 = [-50 + rand(2,1)*100;0];

Pars.fc = 1e9;  %carrier frequency
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;

modOrder=8;
FFTLength=64;
NumSymbols=100;
snr=12;


% definizione MIMO array
Geometry.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');


%signal generation
[ofdmMod,waveform_t1,in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder); %signal tx by UE 1
[~,waveform_t2,in_t2]=OFDMsignal(FFTLength, NumSymbols,modOrder);%signal tx by UE 2
[~,waveform_i,in_i]=OFDMsignal(FFTLength, NumSymbols,modOrder); %signal tx by interf.

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
%showResourceMapping(ofdmMod);

%OFDM demodulators
title='TX Constellation by UE';
bits=OFDMDemod(ofdmMod,waveform_t1,modOrder,true,title);
[numbError,ratio]=biterr(in_t1,bits);

polarLMS=figure('Name','beams');
ueplot=figure('Name','Positions');



%% while cycle
while(1)
    %movements
    t1 = t1+[-10+20*rand(2,1);0]; %random movements (between -5 and 5)
    t2 = t2+[-10+20*rand(2,1);0];
    i1 = i1+[-10+20*rand(2,1);0];
   i2 = i2+[-10+20*rand(2,1);0];
    %position plot
    v=[t1,t2,i1,i2];
    figure(ueplot);
    plot3(bs(1),bs(2),bs(3),'vr','MarkerSize',9,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);
    hold on;
    plot3(v(1,(1)),v(2,(1)),v(3,(1)),'ob','LineWidth',1.5);
    hold on
    plot3(v(1,(2)),v(2,(2)),v(3,(2)),'og','LineWidth',1.5)%tracked terminals
    plot3(v(1,(3:4)),v(2,(3:4)),v(3,(3:4)),'x', 'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); %interferers
    xlim([-100,100])
    ylim([-100,100])
    zlim([0,50])
    grid on;
    legend('Base Station','Terminal 1','Terminal 2','Interferers');
    hold off;    


    % compute azimuth and elevation angles for each terminal (UE+interf.)
    az_t1=rad2deg(atan2(t1(1),t1(2)));
    az_t2=rad2deg(atan2(t2(1),t2(2)));
    az_i1=rad2deg(atan2(i1(1),i1(2)));
    az_i2=rad2deg(atan2(i2(1),i2(2)));


    el_t1=rad2deg(atan2(bs(3),sqrt(t1(1)^2+t1(2)^2)));
    el_t2=rad2deg(atan2(bs(3),sqrt(t2(1)^2+t2(2)^2)));
    el_i1=rad2deg(atan2(bs(3),sqrt(i1(1)^2+i1(2)^2)));
    el_i2=rad2deg(atan2(bs(3),sqrt(i2(1)^2+i2(2)^2)));
    
    %distances from base station in xy plane (BS is in (0,0))
    t1_dist = sqrt(t1(1)^2+t1(2)^2);
    i1_dist = sqrt(i1(1)^2+i1(2)^2);
    t2_dist = sqrt(t2(1)^2+t2(2)^2);
    i2_dist = sqrt(i2(1)^2+i2(2)^2);
    
    %distances travelled by the signal
    t1_dist_BS=sqrt(bs(3)^2+t1_dist^2);
    t2_dist_BS=sqrt(bs(3)^2+t2_dist^2);
    i1_dist_BS=sqrt(bs(3)^2+i1_dist^2);
    i2_dist_BS=sqrt(bs(3)^2+i2_dist^2);
    

    % pathloss
    path_loss_t1 = ((4*pi*t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*t2_dist_BS)/Pars.lambda)^2;
    path_loss_i1 = ((4*pi*i1_dist_BS)/Pars.lambda)^2;
    path_loss_i2 = ((4*pi*i2_dist_BS)/Pars.lambda)^2;
    
 
    t1Angles = [az_t1 el_t1];
    t2Angles = [az_t2 el_t2];
    i1Angles = [az_i1 el_i1];
    i2Angles = [az_i2 el_i2];
    
    %% Channel generation

    receivedW = collectPlaneWave(Geometry.BSarray, [waveform_t1*(1/sqrt(path_loss_t1)) waveform_t2*(1/sqrt(path_loss_t2)) waveform_i*(1/sqrt(path_loss_i1)) waveform_i*(1/sqrt(path_loss_i2))], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);
  
    %add noise
    chOut = awgn(receivedW, snr, 'measured');

    %Compute BER without BF on UE1
    title='Constellation without beamforming on UE1';
    bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors no beam 1: %d',numbError);
    fprintf('\nBER no beam 1: %.2f',ratio);
     
     %Compute BER without BF on UE2
     title='Constellation without beamforming on UE2';
     bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,true,title);
     [numbError,ratio]=biterr(in_t2,bits);
     fprintf('\nNumber of errors no beam 2: %d',numbError);
     fprintf('\nBER no beam 2: %.2f',ratio);
  %% phase shift BF
  
    % beamforming UE1 with PhaseShiftBeamformer (PSB)
    beamformerV1 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t1Angles','WeightsOutputPort',true);
    [y1,w1] = beamformerV1(chOut);

     %Compute BER with PSB BF on UE1
    title='Constellation with PhaseShiftBeamformer on UE1';
    bits=OFDMDemod(ofdmMod,y1,modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors phase shift 1: %d',numbError);
    fprintf('\nBER phase shift 1: %f',ratio);
     
   % beamforming UE2 with PSB
    beamformerV2 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t2Angles','WeightsOutputPort',true);
    [y2,w2] = beamformerV2(chOut);

     %Compute BER with PSB BF on UE2
    title='Constellation with PhaseShiftBeamformer on UE2';
    bits=OFDMDemod(ofdmMod,y2,modOrder,true,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2: %d',numbError);
    fprintf('\nBER phase shift 2: %f',ratio);
    
    %% plots of radiation pattern with PSB

    figure(polarLMS);
    polarplot( deg2rad(az_t1),t1_dist_BS/max(t1_dist_BS,t2_dist_BS), 'or','LineWidth',1.5)
    hold on
    polarplot( deg2rad(az_t2),t2_dist_BS/max(t1_dist_BS,t2_dist_BS),'ob','LineWidth',1.5)
    hold on

    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',conj(w1));
    hold on 
    polarplot(H,'r')

    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',conj(w2));

    hold on 
    polarplot(H,'b')
   
   
    
    %% LMS algorithm
   
    %compute LMS weights 
    optimalWeight1 = LMSalgorithm(chOut,waveform_t1,numArrayElements);  
    optimalWeight2 = LMSalgorithm(chOut,waveform_t2,numArrayElements);   

    
    %multiply rx signal by weiths  
    y1=chOut*((optimalWeight1));
    y2=chOut*((optimalWeight2));     
    
    %Compute BER with LMS BF on UE1
    title='Constellation with LMS on UE1';
     bits=OFDMDemod(ofdmMod,(y1),modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
    fprintf('\nBER LMS1: %f',ratio);
    
    %Compute BER with LMS BF on UE1
    title='Constellation with LMS on UE2';
    bits=OFDMDemod(ofdmMod,(y2),modOrder,true,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors LMS2: %d',numbError);
    fprintf('\nBER LMS2: %f\n',ratio);
    
    %% plots of radiation pattern with LMS BF
    figure(polarLMS);
    H1= pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type','power','CoordinateSystem','polar','Weights',(optimalWeight1));

    hold on
    polarplot(H1,'--r')

    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type','power','CoordinateSystem','polar','Weights',(optimalWeight2)) ;
    
    hold on 
    polarplot(H,'--b')
    legend ('UE1', 'UE2','array directivity pattern with PhaseShift BF on UE1','array directivity pattern with PhaseShift BF on UE2','array directivity pattern with LMS BF on UE1','array directivity pattern with LMS BF on UE2','Location','northeastoutside')
        
    hold off
    pause(1);

end