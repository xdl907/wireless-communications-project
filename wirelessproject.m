%random coordinate generation
bs = [0,0,50]; %base station
t1 = [-50 + rand(2,1)*100;0]; %terminal1
t2 = [-50 + rand(2,1)*100;0]; %terminal2
i1 = [-50 + rand(2,1)*100;0]; %interferers
i2 = [-50 + rand(2,1)*100;0];

while(1)
    t1 = t1+[-5+10*rand(2,1);0]; %random movements (between -5 and 5)
    t2 = t2+[-5+10*rand(2,1);0];
    i1 = i1+[-5+10*rand(2,1);0];
    i2 = i2+[-5+10*rand(2,1);0];
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


    %azimuth and elevation

    az_t1=rad2deg(atan2(t1(1),t1(2)));
    az_t2=rad2deg(atan2(t2(1),t2(2)));
    az_i1=rad2deg(atan2(i1(1),i1(2)));
    az_i2=rad2deg(atan2(i2(1),i2(2)));


    el_t1=rad2deg(atan2(bs(3),sqrt(t1(1)^2+t1(2)^2)));
    el_t2=rad2deg(atan2(bs(3),sqrt(t2(1)^2+t2(2)^2)));
    el_i1=rad2deg(atan2(bs(3),sqrt(i1(1)^2+i1(2)^2)));
    el_i2=rad2deg(atan2(bs(3),sqrt(i2(1)^2+i2(2)^2)));
    
    % Initialize system constants
    gc = helperGetDesignSpecsParameters();
    a = rand(1)*100
    % Tunable parameters
    tp.txPower = 9;           % watt
    tp.txGain = -8;           % dB
    tp.mobileRange = 2750;    % m
    tp.mobileAngle = az_t1;       % degrees
    tp.interfPower = 1;       % watt
    tp.interfGain = -20;      % dB
    tp.interfRange = 9000;    % m
    tp.interfAngle =   az_i1;    % degrees
    tp.numTXElements = 8;       
    tp.steeringAngle = az_t1;     % degrees
    tp.rxGain = 108.8320 - tp.txGain; % dB

    numTx= tp.numTXElements;
    
    %% abbiamo cambiato backbaffled da true a false
    %% Signal Transmission 
    % First, configure the system's transmitter. 

    [encoder,scrambler,modulatorRQAM,modulatorOFDM,steeringvec,transmitter,...
        radiator,pilots,numDataSymbols,frmSz] = helperMIMOTxSetup(gc,tp);

    %%
    % There are many components in the transmitter subsystem, such as the
    % convolutional encoder, the scrambler, the QAM modulator, the OFDM
    % modulator, and so on. The message is first converted to an information
    % bit stream and then passed through source coding and modulation stages to
    % prepare for the radiation.

    txBits = randi([0, 1], frmSz,1);
    coded = encoder(txBits);
    bitsS = scrambler(coded);
    tx = modulatorRQAM(bitsS);

    %%
    % In an OFDM system, the data is carried by multiple sub-carriers that are
    % orthogonal to each other. 

    ofdm1 = reshape(tx, gc.numCarriers,numDataSymbols);

    %%
    % Then, the data stream is duplicated to all radiating elements in the
    % transmitting array

    ofdmData = repmat(ofdm1,[1, 1, numTx]);
    txOFDM = modulatorOFDM(ofdmData, pilots);
    %scale
    txOFDM = txOFDM * ...
        (gc.FFTLength/sqrt(gc.FFTLength-sum(gc.NumGuardBandCarriers)-1));

    % Amplify to achieve peak TX power for each channel
    for n = 1:numTx
        txOFDM(:,n) = transmitter(txOFDM(:,n));
    end

    %%
    % In a MIMO system, it is also possible to separate multiple users spatial
    % division multiplexing (SDMA). In these situations, the data stream is
    % often modulated by a weight corresponding to the desired direction so
    % that once radiated, the signal is maximized in that direction. Because in
    % a MIMO channel, the signal radiated from different elements in an array
    % may go through different propagation environments, the signal radiated
    % from each antenna should be propagated individually. This can be achieved
    % by setting CombineRadiatedSignals to false on the phased.Radiator
    % component.

    radiator.CombineRadiatedSignals = true;

    %%
    % To achieve precoding, the data stream radiated from each antenna in the
    % array is modulated by a phase shift corresponding to its radiating
    % direction. The goal of this precoding is to ensure these data streams add
    % in phase if the array is steered toward that direction. Precoding can be
    % specified as weights used at the radiator.

    wR = steeringvec(gc.fc,[-tp.mobileAngle;0]);

    %% 
    % Meanwhile, the array is also steered toward a given steering angle, so
    % the total weights are a combination of both precoding and the steering
    % weights.

    wT = steeringvec(gc.fc,[tp.steeringAngle;0]);
    weight = wT.* wR;

    %%
    % The transmitted signal is thus given by

    txOFDM = radiator(txOFDM,repmat([tp.mobileAngle;0],1,numTx),conj(weight));

    %%
    % Note that the transmitted signal, txOFDM, is a matrix whose columns
    % represent data streams radiated from the corresponding elements in the
    % transmit array.
    %%
    % The entire scene can be depicted in the figure below
    figure(2);

    
    helperPlotMIMOEnvironment(gc, tp);

    hold off
    pause(1)
end

