clear, close all;
%% general settings 
% Random coordinate generation
g.bs = [0,0,50]; %base station
g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal 1
g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal 2
g.i1 = [-50 + rand(2,1)*100;1.5]; %interferers (static)
g.i2 = [-50 + rand(2,1)*100;1.5];


Pars.fc = 1e9; %carrier frequency (add LTE one?)
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;
Pars.SNR = 20;
numArrayElements=4;

Fsin = 600;
Fsin2 = 500;
Fsin3 = 750;
Ts = 1e-5;
Fsample = 1 / Ts;
TsVect_t1 = 0:Ts:5/Fsin;
TsVect_t2 = 0:Ts:5/Fsin;
TsVect_i = 0:Ts:5/Fsin;

%tx signals
sinusoid_waveform=sin(2*pi*Fsin3*TsVect_i);
sinusoid_waveform(2,:)=sin(2*pi*Fsin3*TsVect_i);
sinusoid_waveform(3,:) = sin(2*pi*Fsin*TsVect_t1);
sinusoid_waveform(4,:)=sin(2*pi*Fsin2*TsVect_t2);
waveform_LEN=size(sinusoid_waveform(1,:));

% MIMO array
g.BSarray = phased.URA('Size', [numArrayElements numArrayElements], 'ElementSpacing', [Pars.lambda/2 Pars.lambda/2], 'ArrayNormal', 'x');
g.bs_antenna_pos=getElementPosition(g.BSarray); %positions of the antennas in the array

%% Array definitions for the Quadriga channel model
l = qd_layout;
l.set_scenario('QuaDRiGa_UD2D_LOS'); %propagation model

array_tx = qd_arrayant('omni');
array_rx = qd_arrayant('omni');
array_rx.no_elements = numArrayElements*numArrayElements;
array_rx.element_position=g.bs_antenna_pos; %fetching the URA antenna positions for the quadriga layout

%% Quadriga channel model
l.tx_array = array_tx; %array definitions
l.rx_array = array_rx;
l.no_rx = 1;
l.no_tx = 4;
while (1)

tx_track1= qd_track('linear',20*rand,rand*2*pi); %track of the 1st terminal
tx_track2= qd_track('linear',30*rand,rand*2*pi); %track of the 2nd terminal 
tx_track1.name = 'UE1';
tx_track2.name = 'UE2';

l.tx_track(1,1) = copy(tx_track1); %track establishment
l.tx_track(1,2) = copy(tx_track2);
l.rx_position = g.bs'; %position fetching
l.tx_position=[g.t1,g.t2,g.i1,g.i2];
l.visualize();
l.set_pairing; %pair all tx/rx links
chan = l.get_channels(); %channel generation

%update current location
g.t1(1)=g.t1(1)+tx_track1.positions(1,2);
g.t1(2)=g.t1(2)+tx_track1.positions(2,2);
g.t1(3)=g.t1(3)+tx_track1.positions(3,2);

g.t2(1)=g.t2(1)+tx_track2.positions(1,2);
g.t2(2)=g.t2(2)+tx_track2.positions(2,2);
g.t2(3)=g.t2(3)+tx_track2.positions(3,2);

%% channels convolution
chTaps=size(chan(3).delay); % channel 1&2 -> interf_Rx    channel 3 -> UE1_Rx    channel 4 -> UE2_Rx    
chOut=zeros(chTaps(1),waveform_LEN(2),4);
TsVect=0:Ts:Ts*(waveform_LEN(2)-1);

for channel=1:1:4
    for antenna=1:1:chTaps(1)
        for path=1:1:chTaps(3)

            inX=TsVect-chan(channel).delay(antenna,1,path,1);
            inY=interp1(TsVect,sinusoid_waveform(channel,:),inX,'pship');
            chOut(antenna,:,channel)=inY*chan(channel).coeff(antenna,1,path,1)+chOut(antenna,:,channel);

        end
    end
end
%sum responses of all single channels in one
chOut(:,:,1)=chOut(:,:,1)+chOut(:,:,2)+chOut(:,:,3)+chOut(:,:,4);
chOut(:,:,2:4)=[];

figure
hold on
plot(sinusoid_waveform(3,:))
for i=1:3
    chOut(i,:)=chOut(i,:)/max((chOut(i,:)));%normalize values
    plot (real(chOut(i,:))) 
end
legend ('tx sinusoid', 'rx signal antenna element 1', 'rx signal antenna element 2', 'rx signal antenna element 3');
hold off

chOut = awgn(chOut, Pars.SNR, 'measured');

pause(2)

end







