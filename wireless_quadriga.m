%random coordinate generation
g.bs = [0,0,25]; %base station
g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal1
g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal2
g.i1 = [-50 + rand(2,1)*100;1.5]; %interferers
g.i2 = [-50 + rand(2,1)*100;1.5];

%phased array definition
Pars.fc = 1e9;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;
g.BSarray = phased.URA('Size',[4 4], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');
g.bs_antenna_pos=getElementPosition(g.BSarray);

%array definitions for the quadriga channel model
l = qd_layout;
l.set_scenario('QuaDRiGa_UD2D_LOS');

array_tx = qd_arrayant('omni');
array_rx = qd_arrayant('omni');
array_rx.no_elements = 16;
array_rx.element_position=g.bs_antenna_pos;

%Quadriga channel model
l.tx_array = array_tx;
l.rx_array = array_rx;
l.no_rx = 1;
l.no_tx = 4;
tx_track1= qd_track('linear',30,pi);
tx_track2= qd_track('linear',20,pi/2);
tx_track1.name = 'UE1';
tx_track2.name = 'UE2';

l.tx_track(1,1) = copy(tx_track1);
l.tx_track(1,2) = copy(tx_track2);
l_rx_position = g.bs';
l.tx_position=[g.t1,g.t2,g.i1,g.i2];
l.visualize();
