clear all;
clc;
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

wStart=complex(ones(numArrayElements*numArrayElements,1));

% definizione MIMO array
Geometry.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');
% Geometry.BSarray.Element.FrequencyRange = [50.0e6 1000.0e6];

%generazione dei segnali
Fsin = 600;
Fsin2 = 500;
% Fsin3 = 750;
Ts = 1e-5;
Fsample = 1 / Ts;
TsVect_t1 = 0:Ts:5/Fsin;
TsVect_t2 = 0:Ts:5/Fsin;
% TsVect_i = 0:Ts:5/Fsin;
sinusoid_waveform_t1 = sin(2*pi*Fsin*TsVect_t1);
sinusoid_waveform_t2 = sin(2*pi*Fsin2*TsVect_t2);
% sinusoid_waveform_i = sin(2*pi*Fsin3*TsVect_i);

%% calcolo dei dati utili
%plot procedure
    v=[t1,t2,i1,i2];
    figure(1)
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
    
%% generazione del segnale
%     receivedW = collectPlaneWave(Geometry.BSarray, [(sinusoid_waveform_t1*sqrt(path_loss_t1))' (sinusoid_waveform_t2*sqrt(path_loss_t2))' (sinusoid_waveform_i*sqrt(path_loss_i1))' (sinusoid_waveform_i*sqrt(path_loss_i2))'], [t1Angles' t2Angles' i1Angles' i2Angles'], Pars.fc);  %% dobbiamo aggiungere l'interferenza e moltiplicare per la radice quadrata del pathloss
    receivedW = collectPlaneWave(Geometry.BSarray, [(sinusoid_waveform_t1)' (sinusoid_waveform_t2)'], [t1Angles' t2Angles'], Pars.fc);    
    Pars.SNR = 20;
    chOut = awgn(receivedW, Pars.SNR, 'measured'); %% segnale del segnale in entrata alla BS

%% generazione DOA algo

estimator = phased.MUSICEstimator2D('SensorArray', Geometry.BSarray,...
    'OperatingFrequency', Pars.fc, 'ForwardBackwardAveraging', true, 'NumSignalsSource', 'Property',...
    'DOAOutputPort', true, 'NumSignals', 2, 'AzimuthScanAngles', -90:0.5:90, ...
    'ElevationScanAngles', -90:0.5:90); % problema della generazione degli angoli: oltre l'intervallo +-90
                                        % diventa incasinato gestirselo

[~,doas] = estimator(chOut);

disp(doas);

%plot lo spazio in 2D
figure
plotSpectrum(estimator);


%% beamforming con i valori angolari effettivi
beamformer1s = phased.PhaseShiftBeamformer('SensorArray', Geometry.BSarray, 'OperatingFrequency', Pars.fc,...
    'Direction', t1Angles, 'WeightsOutputPort', true);
[arrOut1, w1] = beamformer1s(chOut);

beamformer2s = phased.PhaseShiftBeamformer('SensorArray', Geometry.BSarray, 'OperatingFrequency', Pars.fc,...
    'Direction', t2Angles, 'WeightsOutputPort', true);
[arrOut2, w2] = beamformer2s(chOut);

%% beamforming con i valori angolari stimati con DOA
doas1 = doas(:,2);
doas2 = doas(:,1);

beamformer1e = phased.PhaseShiftBeamformer('SensorArray', Geometry.BSarray, 'OperatingFrequency', Pars.fc,...
    'Direction', doas1, 'WeightsOutputPort', true);
[arrOut1e, w1e] = beamformer1e(chOut);

beamformer2e = phased.PhaseShiftBeamformer('SensorArray', Geometry.BSarray, 'OperatingFrequency', Pars.fc,...
    'Direction', doas2, 'WeightsOutputPort', true);
[arrOut2e, w2e] = beamformer2e(chOut);


%% stampo i grafici polari

 figure(6)
    pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'powerdb','CoordinateSystem','polar','Weights',w1)
      
    hold on 
    
    pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type',...
    'powerdb','CoordinateSystem','polar','Weights',w2)
    hold off;
    
    hold on
    
    pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t1,'PropagationSpeed',Pars.c,'Type',...
    'powerdb','CoordinateSystem','polar','Weights',w1e)
      
    hold on 
    
    pattern(Geometry.BSarray,Pars.fc,[-180:180],el_t2,'PropagationSpeed',Pars.c,'Type',...
    'powerdb','CoordinateSystem','polar','Weights',w2e)
    hold off;
    
    p = polarpattern('gco');
    p.LegendLabels = {'terminal 1','terminal 2'};
