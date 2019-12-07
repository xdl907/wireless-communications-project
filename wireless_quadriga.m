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

%% OFDM signal generation (LTE-A)

%QAM signal generation
M = 64; %modulation order
k = log2(M); %bits per symbol
tx_symbols = randi([0 M-1],2000,1);
refc = qammod(0:M-1,M);
mod_qam = qammod(tx_symbols, M, 'gray');
%add encoding in the future? ask reggiani
constdiag = comm.ConstellationDiagram('ReferenceConstellation',refc, ...
    'XLimits',[-9 9],'YLimits',[-9 9]);
constdiag(mod_qam);
rcv = awgn(mod_qam,10); %AWGN channel for testing
constdiag(rcv);

%OFDM parameters
n_subc = 2048; %subcarriers
cp_len = n_subc*0.25; %with long cyclic prefix in LTE, the length is 25% of the total symbol rate (??? da verificare)





