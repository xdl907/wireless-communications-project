clearvars
close all
rng shuffle
%% environment initialization
%random coordinate generation
bs = [0,0,50]; %base station
t1 = [-50 + rand*100;abs(-50 + rand*100);0]; %terminal1
t2 = [-50 + rand*100;abs(-50 + rand*100);0]; %terminal2
i1 = [-50 + rand*100;abs(-50 + rand*100);0]; %interferers
i2 = [-50 + rand*100;abs(-50 + rand*100);0];


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

% graph initialization
polarLMS=figure('Name','beams');
ueplot=figure('Name','Positions');
aoaSpectrum=figure('Name','AoA spatial spectrum');

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



%% cycle
for c = 1:60
    
    % movements       
    t1(1) = t1(1)-20*rand; 
    t1(2) = t1(2)+20*rand;
    
    t2(1) = t2(1)-20*rand; 
    t2(2) = t2(2)+20*rand;
    
    i1(1) = i1(1)-20*rand; 
    i1(2) = i1(2)+20*rand;
    
    i2(1) = i2(1)-20*rand; 
    i2(2) = i2(2)+20*rand;

    % environment plot generation
    v=[t1,t2,i1,i2];
    figure(ueplot);
     plot3(bs(1),bs(2),bs(3),'vr','MarkerSize',9,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);
    hold on;
    plot3(v(1,(1)),v(2,(1)),v(3,(1)),'ob','LineWidth',1.5);
    hold on
    plot3(v(1,(2)),v(2,(2)),v(3,(2)),'og','LineWidth',1.5)%tracked terminals
    plot3(v(1,(3:4)),v(2,(3:4)),v(3,(3:4)),'x', 'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); %interferers
    zlim([0,50])
    grid on;
    legend('Base Station','Terminal 1','Terminal 2','Interferers');
    hold off;   


    % azimuth and elevation calculation
    az_t1=rad2deg(atan2(t1(1),t1(2)));
    az_t2=rad2deg(atan2(t2(1),t2(2)));
    az_i1=rad2deg(atan2(i1(1),i1(2)));
    az_i2=rad2deg(atan2(i2(1),i2(2)));

    el_t1=rad2deg(atan2(bs(3),sqrt(t1(1)^2+t1(2)^2)));
    el_t2=rad2deg(atan2(bs(3),sqrt(t2(1)^2+t2(2)^2)));
    el_i1=rad2deg(atan2(bs(3),sqrt(i1(1)^2+i1(2)^2)));
    el_i2=rad2deg(atan2(bs(3),sqrt(i2(1)^2+i2(2)^2)));
    
    t1_dist = sqrt(t1(1)^2+t1(2)^2);
    i1_dist = sqrt(i1(1)^2+i1(2)^2);
    t2_dist = sqrt(t2(1)^2+t2(2)^2);
    i2_dist = sqrt(i2(1)^2+i2(2)^2);
    
    t1_dist_BS=sqrt(bs(3)^2+t1_dist^2);
    t2_dist_BS=sqrt(bs(3)^2+t2_dist^2);
    i1_dist_BS=sqrt(bs(3)^2+i1_dist^2);
    i2_dist_BS=sqrt(bs(3)^2+i2_dist^2);
    

    % pathloss calculation
    path_loss_t1 = ((4*pi*t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*t2_dist_BS)/Pars.lambda)^2;
    path_loss_i1 = ((4*pi*i1_dist_BS)/Pars.lambda)^2;
    path_loss_i2 = ((4*pi*i2_dist_BS)/Pars.lambda)^2;
    
    % set of angles
    t1Angles = [az_t1 el_t1];
    t2Angles = [az_t2 el_t2];
    i1Angles = [az_i1 el_i1];
    i2Angles = [az_i2 el_i2];
    
    %% PhaseShiftBeamformer
    % plane wave generation with AWGN and fixed SNR
    receivedW = collectPlaneWave(Geometry.BSarray, [waveform_t1*(1/sqrt(path_loss_t1)) waveform_t2*(1/sqrt(path_loss_t2)) waveform_i*(1/sqrt(path_loss_i1)) waveform_i*(1/sqrt(path_loss_i2))], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);
    Pars.SNR = snr;
    chOut = awgn(receivedW, Pars.SNR, 'measured');
    
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
    'DOAOutputPort', true, 'NumSignals', 4, 'AzimuthScanAngles', -90:0.5:90, ...
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
    fprintf('\nNumber of errors phase shift 1 using DoA: %d',numbError);
    fprintf('\nBER phase shift 1: %f',ratio);
    
    title='Constellation with PhaseShiftBeamformer on UE1 using real angles';
    bits=OFDMDemod(ofdmMod,y1real,modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors phase shift 1 using real angles: %d',numbError);
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
    
    title='Constellation with PhaseShiftBeamformer on UE2 using DoA';
    bits=OFDMDemod(ofdmMod,y2,modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2 for DoA: %d',numbError);
    fprintf('\nBER phase shift 2: %f',ratio);
    
    title='Constellation with PhaseShiftBeamformer on UE2 using real angles';
    bits=OFDMDemod(ofdmMod,y2real,modOrder,false,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2 for real angles: %d',numbError);
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
   
    
    %% LMS algorithm 
   
    optimalWeight1 = LMSalgorithm(chOut,waveform_t1,numArrayElements,Pars.lambda);  
    optimalWeight2 = LMSalgorithm(chOut,waveform_t2,numArrayElements,Pars.lambda);   

    
    for i=1:length(waveform_t1)
        y1(i)=chOut(i,:)*(optimalWeight1);
        y2(i)=chOut(i,:)*(optimalWeight2);
    end
%     
    title='Constellation with LMS on UE1 based on DoA';
     bits=OFDMDemod(ofdmMod,(y1),modOrder,false,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
        fprintf('\nBER LMS1: %f',ratio);
    
    
    title='Constellation with LMS on UE2 based on DoA';
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