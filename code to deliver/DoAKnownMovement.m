clearvars
close all

%% environment initialization
%random coordinate generation
bs = [0,0,50]; %base station
t1 = [abs(-50 + rand(2,1)*100);0]; %terminal1
t2 = [abs(-50 + rand(2,1)*100);0]; %terminal2


Pars.fc = 1e9;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

% MIMO array's parameters generation
numArrayElements=4;

modOrder=8;
FFTLength=64;
NumSymbols=1000;
snr=15;


% MIMO array definition
Geometry.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');

% signal generation

[ofdmMod,waveform_t1,in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder);
[~,waveform_t2,in_t2]=OFDMsignal(FFTLength, NumSymbols,modOrder);
[~,waveform_i,in_i]=OFDMsignal(FFTLength, NumSymbols,modOrder);

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
% showResourceMapping(ofdmMod_t1);

%OFDM demodulators
title='TX Constellation by UE';
bits=OFDMDemod(ofdmMod,waveform_t1,modOrder,false,title);
% [numbError,ratio]=biterr(in_t1,bits);
% numbError
% ratio


% graph initialization
polarLMS=figure('Name','beams');
ueplot=figure('Name','Positions');
aoaSpectrum=figure('Name','AoA spatial spectrum');

%step for the movement
up1 = 0;
up2 = 0;


%% cycle
for c = 1:60
    % movements
    if (t1(1)>= 100 || t1(2) >= 100)
        up1 = 20;
    end
    if (up1>0)
        t1 = abs(t1+[-2;-1;0]);
        up1 = up1 -1;
    else
        t1 = abs(t1+[2;1;0]);
    end
    if (t2(1) >= 100 || t2(2) >= 100)
        up2 = 20;
    end
    if (up2>0)
        t2 = abs(t2+[-2;-1;0]);
        up2 = up2 - 1;
    else
        t2 = abs(t2+[2;1;0]);
    end

    
    % plot generation
    v=[t1,t2];
    figure(ueplot);
  plot3(bs(1),bs(2),bs(3),'vr','MarkerSize',9,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);
    hold on;
    plot3(v(1,(1)),v(2,(1)),v(3,(1)),'ob','LineWidth',1.5);
    hold on
    plot3(v(1,(2)),v(2,(2)),v(3,(2)),'og','LineWidth',1.5)%tracked terminals
    xlim([-100,100])
    ylim([-100,100])
    zlim([0,50])
    grid on;
    legend('Base Station','Terminal 1','Terminal 2');
    hold off;


    % azimuth and elevation calculation

    az_t1=rad2deg(atan2(t1(1),t1(2)));
    az_t2=rad2deg(atan2(t2(1),t2(2)));


    el_t1=rad2deg(atan2(bs(3),sqrt(t1(1)^2+t1(2)^2)));
    el_t2=rad2deg(atan2(bs(3),sqrt(t2(1)^2+t2(2)^2)));

    
    t1_dist_BS = sqrt(t1(1)^2+t1(2)^2);
    t2_dist_BS = sqrt(t2(1)^2+t2(2)^2);
   
    

   % pathloss calculation
    path_loss_t1 = ((4*pi*t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*t2_dist_BS)/Pars.lambda)^2;

    
    % set of angles
    t1Angles = [az_t1 el_t1];
    t2Angles = [az_t2 el_t2];
    
    %% PhaseShiftBeamformer
    % plane wave generation with AWGN and fixed SNR
    receivedW = collectPlaneWave(Geometry.BSarray, [waveform_t1*(1/sqrt(path_loss_t1)) waveform_t2*(1/sqrt(path_loss_t2))], [t1Angles' t2Angles'], Pars.fc);
    Pars.SNR = snr;
    chOut = awgn(receivedW, Pars.SNR, 'measured'); %% segnale del segnale in entrata alla BS
    
    % constellation without beamforming
    title='Constellation without beamforming on UE1';
    bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors no beam 1: %d',numbError);
        fprintf('\nBER no beam 1: %.2f',ratio);
    
     title='Constellation without beamforming on UE2';
     bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
      fprintf('\nNumber of errors no beam 2: %d',numbError);
        fprintf('\nBER no beam 2: %.2f',ratio);
   
    %% Direction-of-Arrival algorithm
    % estimator generation
    estimator = phased.MUSICEstimator2D('SensorArray', Geometry.BSarray,...
    'OperatingFrequency', Pars.fc, 'ForwardBackwardAveraging', true, 'NumSignalsSource', 'Property',...
    'DOAOutputPort', true, 'NumSignals', 2, 'AzimuthScanAngles', -90:0.5:90, ...
    'ElevationScanAngles', -90:0.5:90);

    % estimation of DoA
    [~,doas] = estimator(chOut);
    
    % angles assignment
    doa1=doas(:,1);
    doa2=doas(:,2);
    
    % spectrum plot
    figure(aoaSpectrum)
    plotSpectrum(estimator);

    %% Beamforming and performance metrics
    % beamforming  of the first UE with DoA angles
    beamformerV1 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',doa1,'WeightsOutputPort',true);
    [y1,w1] = beamformerV1(chOut);
    
    % beamforming  of the first UE with real angles
    beamformerV1real = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t1Angles','WeightsOutputPort',true);
    [y1real,w1real] = beamformerV1real(chOut);
    
    % related costellations
    title='Constellation with PhaseShiftBeamformer on UE1 using DoA';
    bits=OFDMDemod(ofdmMod,y1,modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors phase shift 1: %d',numbError);
    fprintf('\nBER phase shift 1 DoA: %f',ratio);
    
    title='Constellation with PhaseShiftBeamformer on UE1 using real angles';
    bits=OFDMDemod(ofdmMod,y1real,modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors phase shift 1: %d',numbError);
    fprintf('\nBER phase shift 1: %f',ratio);
 
    % beamforming  of the second UE with DoA angles
    beamformerV2 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',doa2,'WeightsOutputPort',true);
    [y2,w2] = beamformerV2(chOut);
    
    % beamforming  of the second UE with real angles
    beamformerV2real = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t2Angles','WeightsOutputPort',true);
    [y2real,w2real] = beamformerV2real(chOut);
    
    % related costellations
    title='Constellation with PhaseShiftBeamformer on UE2 using DoA';
    bits=OFDMDemod(ofdmMod,y2,modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2: %d',numbError);
    fprintf('\nBER phase shift 2 DoA: %f',ratio);
    
    title='Constellation with PhaseShiftBeamformer on UE2 using real angles';
    bits=OFDMDemod(ofdmMod,y2real,modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2: %d',numbError);
    fprintf('\nBER phase shift 2: %f',ratio);
    
    % plot of the beam's graphs
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
    
    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',conj(w1real));
    hold on 
    polarplot(H,':r')
    
    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',conj(w2real));
    hold on 
    polarplot(H,':b')

    
    %% algoritmo LMS 
   
    optimalWeight1 = LMSalgorithm(chOut,waveform_t1,numArrayElements,Pars.lambda);  
    optimalWeight2 = LMSalgorithm(chOut,waveform_t2,numArrayElements,Pars.lambda);   

    
    for i=1:length(waveform_t1)
        y1(i)=chOut(i,:)*(optimalWeight1);
        y2(i)=chOut(i,:)*(optimalWeight2);
    end
%     
    title='Constellation with LMS on UE1';
     bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
        fprintf('\nBER LMS1: %f',ratio);
    
    
    title='Constellation with LMS on UE2';
     bits=OFDMDemod(ofdmMod,(y2),modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
       fprintf('\nNumber of errors LMS2: %d',numbError);
        fprintf('\nBER LMS2: %f\n',ratio);
    

    % polar plot of the results of LMS algorithm.
    figure(polarLMS);
   H1= pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',(optimalWeight1));

   hold on
    polarplot(H1,'--r')

    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',(optimalWeight2)) ;

    hold on 
    polarplot(H,'--b')
    legend ('UE1', 'UE2','array directivity pattern with PhaseShift BF on UE1 doa1','array directivity pattern with PhaseShift BF on UE2 doa2','array directivity pattern with PhaseShift BF on UE1 based on real angles', 'array directivity pattern with PhaseShift BF on UE2 based on real angles','array directivity pattern with LMS BF on UE1','array directivity pattern with LMS BF on UE2','Location','northeastoutside')

    hold off
    pause(2);

end