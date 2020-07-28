clearvars
close all

%% impostazione dell'ambiente
%random coordinate generation
bs = [0,0,50]; %base station
t1 = [-50 + rand(2,1)*100;0]; %terminal1
t2 = [-50 + rand(2,1)*100;0]; %terminal2
i1 = [-50 + rand(2,1)*100;0]; %interferers
i2 = [-50 + rand(2,1)*100;0];

Pars.fc = 1e9;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;

modOrder=8;
FFTLength=64;
NumSymbols=1000;
snr=15;


% definizione MIMO array
Geometry.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');

%generazione dei segnali

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
bits=OFDMDemod(ofdmMod,waveform_t1,modOrder,true,title);
% [numbError,ratio]=biterr(in_t1,bits);
% numbError
% ratio


    
% polarPhased=figure('Name','beam phased shift');
polarLMS=figure('Name','beams');
ueplot=figure('Name','Positions');



%% inizio del ciclo while
while(1)
    % spostamenti
    t1 = t1+[-10+20*rand(2,1);0]; %random movements (between -5 and 5)
    t2 = t2+[-10+20*rand(2,1);0];
    i1 = i1+[-10+20*rand(2,1);0];
    i2 = i2+[-10+20*rand(2,1);0];
    %plot procedure
    v=[t1,t2,i1,i2];
    figure(ueplot);
    plot3(bs(1),bs(2),bs(3),'vr','MarkerSize',9,'LineWidth',1.5,'Color',[0.6350 0.0780 0.1840]);
    hold on;
    plot3(v(1,(1:2)),v(2,(1:2)),v(3,(1:2)),'ob','LineWidth',1.5); %tracked terminals
    plot3(v(1,(3:4)),v(2,(3:4)),v(3,(3:4)),'x', 'Color',[0.8500 0.3250 0.0980],'LineWidth',1.5); %interferers
    xlim([-100,100])
    ylim([-100,100])
    zlim([0,50])
    grid on;
    legend('Base Station','Terminals','Interferers');
    hold off;    


    % calcolo azimuth and elevation

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
    

    % calcolo del pathloss
    path_loss_t1 = ((4*pi*t1_dist_BS)/Pars.lambda)^2;
    path_loss_t2 = ((4*pi*t2_dist_BS)/Pars.lambda)^2;
    path_loss_i1 = ((4*pi*i1_dist_BS)/Pars.lambda)^2;
    path_loss_i2 = ((4*pi*i2_dist_BS)/Pars.lambda)^2;
    
 
    t1Angles = [az_t1 el_t1];
    t2Angles = [az_t2 el_t2];
    i1Angles = [az_i1 el_i1];
    i2Angles = [az_i2 el_i2];
    
    %% PhaseShiftBeamformer

%     receivedW = collectPlaneWave(Geometry.BSarray, [(sinusoid_waveform_t1/sqrt(path_loss_t1))' (sinusoid_waveform_t2/sqrt(path_loss_t2))' (sinusoid_waveform_i/sqrt(path_loss_i1))' (sinusoid_waveform_i/sqrt(path_loss_i2))'], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);  %% dobbiamo aggiungere l'interferenza e moltiplicare per la radice quadrata del pathloss
    receivedW = collectPlaneWave(Geometry.BSarray, [waveform_t1*(1/sqrt(path_loss_t1)) waveform_t2*(1/sqrt(path_loss_t2)) waveform_i*(1/sqrt(path_loss_i1)) waveform_i*(1/sqrt(path_loss_i2))], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);
    Pars.SNR = snr;
    chOut = awgn(receivedW, Pars.SNR, 'measured'); %% segnale del segnale in entrata alla BS
    
    title='Constellation without beamforming on UE1';
    bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors no beam 1: %d',numbError);
        fprintf('\nBER no beam 1: %.2f',ratio);
    
     title='Constellation without beamforming on UE1';
     bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,true,title);
    [numbError,ratio]=biterr(in_t2,bits);
      fprintf('\nNumber of errors no beam 2: %d',numbError);
        fprintf('\nBER no beam 2: %.2f',ratio);
   
    % beamforming UE1 
    beamformerV1 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t1Angles','WeightsOutputPort',true);
    [y1,w1] = beamformerV1(chOut);
    
    title='Constellation with PhaseShiftBeamformer on UE1';
    bits=OFDMDemod(ofdmMod,y1,modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors phase shift 1: %d',numbError);
    fprintf('\nBER phase shift 1: %f',ratio);
 
    % beamforming UE2
    beamformerV2 = phased.PhaseShiftBeamformer('SensorArray',Geometry.BSarray,...
    'OperatingFrequency',Pars.fc,'PropagationSpeed',Pars.c,...
    'Direction',t2Angles','WeightsOutputPort',true);
    [y2,w2] = beamformerV2(chOut);
    
    title='Constellation with PhaseShiftBeamformer on UE2';
    bits=OFDMDemod(ofdmMod,y2,modOrder,true,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors phase shift 2: %d',numbError);
    fprintf('\nBER phase shift 2: %f',ratio);
    
    
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
   
   
    
    %% algoritmo LMS 
   
    optimalWeight1 = LMSalgorithm(chOut,waveform_t1,numArrayElements,Pars.lambda);  
    optimalWeight2 = LMSalgorithm(chOut,waveform_t2,numArrayElements,Pars.lambda);   

    
    for i=1:length(waveform_t1)
        y1(i)=chOut(i,:)*(optimalWeight1);
        y2(i)=chOut(i,:)*(optimalWeight2);
    end
    
    title='Constellation with LMS on UE1';
     bits=OFDMDemod(ofdmMod,(y1),modOrder,true,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
        fprintf('\nBER LMS1: %f',ratio);
    
    
    title='Constellation with LMS on UE2';
     bits=OFDMDemod(ofdmMod,(y2),modOrder,true,title);
    [numbError,ratio]=biterr(in_t2,bits);
       fprintf('\nNumber of errors LMS2: %d',numbError);
        fprintf('\nBER LMS2: %f\n',ratio);
    
    % plot dei risultati.
    figure(polarLMS);
   H1= pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',(optimalWeight1));

   hold on
    polarplot(H1,'--r')

    H=pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type',...
    'power','CoordinateSystem','polar','Weights',(optimalWeight2)) ;

    hold on 
    polarplot(H,'--b')
    legend ('UE1', 'UE2','array directivity pattern with PhaseShift BF on UE1','array directivity pattern with PhaseShift BF on UE2','array directivity pattern with LMS BF on UE1','array directivity pattern with LMS BF on UE2','Location','northeastoutside')
        
    hold off
    pause(5);

end