clear, close all;
%% Random coordinate generation
g.bs = [0,0,25]; %base station
g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal 1
g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal 2
g.i1 = [-50 + rand(2,1)*100;1.5]; %interferers
g.i2 = [-50 + rand(2,1)*100;1.5];

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

%% Generate Uplink RMC Waveform
% Uplink RMC configuration:
cfg = struct('RC', 'A5-7', ...
    'NULRB', 100, ...
    'DuplexMode', 'FDD', ...
    'NCellID', 0, ...
    'RNTI', 1, ...
    'TotSubframes', 10, ...
    'Windowing', 0);
cfg.PUSCH.RVSeq = [0 2 3 1];
antennaPort = 1;
cfg = lteRMCUL(cfg);
% input bit source:
in = [1; 0; 0; 1];
% waveform generation:
[waveform, grid, cfg] = lteRMCULTool(cfg, in);
Fs = cfg.SamplingRate; % sample rate of waveform

%% Visualize Uplink RMC Waveform
% Spectrum Analyzer
spectrum = dsp.SpectrumAnalyzer('SampleRate', Fs);
spectrum(waveform);
release(spectrum);




