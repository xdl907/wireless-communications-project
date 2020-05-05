clear, close all;

%% Random coordinate generation
g.bs = [0,0,25]; %base station
g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal 1
g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal 2
g.i1 = [-50 + rand(2,1)*100;1.5]; %interferers
g.i2 = [-50 + rand(2,1)*100;1.5];

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
g.i1_dist = sqrt(g.i1(1)^2+g.i1(2)^2);
g.t2_dist = sqrt(g.t2(1)^2+g.t2(2)^2);
g.i2_dist = sqrt(g.i2(1)^2+g.i2(2)^2);

%signal path
g.t1_dist_BS=sqrt(g.bs(3)^2+g.t1_dist^2);
g.t2_dist_BS=sqrt(g.bs(3)^2+g.t2_dist^2);
g.i1_dist_BS=sqrt(g.bs(3)^2+g.i1_dist^2);
g.i2_dist_BS=sqrt(g.bs(3)^2+g.i2_dist^2);
    
%% Phased array definition
Pars.fc = 1e9; %carrier frequency (add LTE one?)
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc; %wavelength 
g.BSarray = phased.URA('Size',[4 4], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x'); %URA definition
g.bs_antenna_pos=getElementPosition(g.BSarray); %positions of the antennas in the array

%% Array definitions for the Quadriga channel model
l = qd_layout;
l.set_scenario('QuaDRiGa_UD2D_LOS'); %propagation model

array_tx = qd_arrayant('omni');
array_rx = qd_arrayant('omni');
array_rx.no_elements = 16;
array_rx.element_position=g.bs_antenna_pos; %fetching the URA antenna positions for the quadriga layout

%% Quadriga channel model
l.tx_array = array_tx; %array definitions
l.rx_array = array_rx;
l.no_rx = 1;
l.no_tx = 4;
tx_track1= qd_track('linear',30,pi); %track of the 1st terminal
tx_track2= qd_track('linear',20,pi/2); %track of the 2nd terminal (add random values?)
tx_track1.name = 'UE1';
tx_track2.name = 'UE2';

l.tx_track(1,1) = copy(tx_track1); %track establishment
l.tx_track(1,2) = copy(tx_track2);
l.rx_position = g.bs'; %position fetching
l.tx_position=[g.t1,g.t2,g.i1,g.i2];
l.visualize();
l.set_pairing; %pair all tx/rx links
chan = l.get_channels(); %channel generation

%to be continued...

%% OFDM signal generation 
load('ofdmwaveform.mat');
% Spectrum Analyzer
spectrum = dsp.SpectrumAnalyzer('SampleRate', waveStruct.Fs);
spectrum(waveStruct.waveform);
release(spectrum);

% OFDM Subcarrier Mapping
ofdmMod = comm.OFDMModulator('FFTLength', 64, ...
    'NumGuardBandCarriers', [6;6], ...
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [16], ...
    'Windowing', false, ...
    'NumSymbols', 100, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [12;26;40;54]);
showResourceMapping(ofdmMod);
ofdmDemod = comm.OFDMDemodulator(ofdmMod);
[dataOut,pilotOut] = step(ofdmDemod,waveStruct.waveform);
constdiag = comm.ConstellationDiagram();

%constdiag(dataOut(:));
constdiag(pilotOut(:));

%% Waveform reception
%pathloss calculation for each terminal
path_loss_t1 = ((4*pi*g.t1_dist_BS)/Pars.lambda)^2;
path_loss_t2 = ((4*pi*g.t2_dist_BS)/Pars.lambda)^2;
path_loss_i1 = ((4*pi*g.i1_dist_BS)/Pars.lambda)^2;
path_loss_i2 = ((4*pi*g.i2_dist_BS)/Pars.lambda)^2;

receivedW = collectPlaneWave(g.BSarray, [(waveStruct.waveform(:)'*sqrt(path_loss_t1))'], ...
    [t1Angles'], Pars.fc);

