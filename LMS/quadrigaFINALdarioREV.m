clear, close all;
%% general settings 
% Random coordinate generation
g.bs = [0,0,50]; %base station
g.t1 = [-50 + rand(2,1)*100;1.5]; %terminal 1
g.t2 = [-50 + rand(2,1)*100;1.5]; %terminal 2
g.i1 = [-50 + rand(2,1)*100;1.5]; %interferers 
g.i2 = [-50 + rand(2,1)*100;1.5];


Pars.fc = 1e9; %carrier frequency
Ts=1/2*Pars.fc;
Pars.c = physconst('LightSpeed');
Pars.lambda = Pars.c/Pars.fc;

numArrayElements=4;

modOrder=8;
FFTLength=64;
NumSymbols=100;
%snr = 10;
graph = false; 
ber1_nobf = [];
ber1_lms = [];
ber2_nobf = [];
ber2_lms = [];
snr_vec = [];
%OFDM tx signals
[~,waveform(:,1),in_i]=OFDMsignal(FFTLength, NumSymbols,modOrder);
waveform(:,2)=waveform(:,1);
[ofdmMod,waveform(:,3),in_t1]=OFDMsignal(FFTLength, NumSymbols,modOrder);
[~,waveform(:,4),in_t2]=OFDMsignal(FFTLength, NumSymbols,modOrder);

% waveform=[waveform_i,waveform_i,waveform_t1,waveform_t2];

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
showResourceMapping(ofdmMod);

%OFDM demodulators
title='TX Constellation by terminals';
bits=OFDMDemod(ofdmMod,waveform(:,3),modOrder,graph,title);
% [numbError,ratio]=biterr(in_t1,bits);
% numbError
% ratio

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
l.tx_array = array_tx; %array definitions
l.rx_array = array_rx;
l.no_rx = 1;
l.no_tx = 4;
l.rx_track.name= 'BS';
l.rx_position = g.bs'; %position fetching
%% Quadriga channel model
%while (1)

tx_track1= qd_track('linear',20*rand,rand*2*pi); %track of the 1st terminal with random length between 0 and 20 m and random direction in [0,2*pi]
tx_track2= qd_track('linear',30*rand,rand*2*pi); %track of the 2nd terminal with random length between 0 and 30 m and random direction in [0,2*pi]
tx_track3= qd_track('linear',10*rand,rand*2*pi); %track of the 1st interf with random length between 0 and 10 m and random direction in [0,2*pi]
tx_track4= qd_track('linear',15*rand,rand*2*pi); %track of the 2nd interf with random length between 0 and 15 m and random direction in [0,2*pi]

tx_track1.name = 'UE1';
tx_track2.name = 'UE2';
tx_track3.name = 'INTERF 1';
tx_track4.name = 'INTERF 2';

l.tx_track(1,1) = copy(tx_track1); %track establishment
l.tx_track(1,2) = copy(tx_track2);
l.tx_track(1,3) = copy(tx_track3); 
l.tx_track(1,4) = copy(tx_track4);

l.tx_position=[g.t1,g.t2,g.i1,g.i2];
pos=l.visualize();
pos.Name='Positions';
pos.NumberTitle='off';
l.set_pairing; %pair all tx/rx links
chan = l.get_channels(); %channel generation

%update current locations
g.t1(1)=g.t1(1)+tx_track1.positions(1,2);
g.t1(2)=g.t1(2)+tx_track1.positions(2,2);
g.t1(3)=g.t1(3)+tx_track1.positions(3,2);

g.t2(1)=g.t2(1)+tx_track2.positions(1,2);
g.t2(2)=g.t2(2)+tx_track2.positions(2,2);
g.t2(3)=g.t2(3)+tx_track2.positions(3,2);

g.i1(1)=g.i1(1)+tx_track3.positions(1,2);
g.i1(2)=g.i1(2)+tx_track3.positions(2,2);
g.i1(3)=g.i1(3)+tx_track3.positions(3,2);

g.i2(1)=g.i2(1)+tx_track4.positions(1,2);
g.i2(2)=g.i2(2)+tx_track4.positions(2,2);
g.i2(3)=g.i2(3)+tx_track4.positions(3,2);

%% channels convolution
for snr=5:15
    chTaps=size(chan(3).delay); % channel 1&2 -> interf_Bs    channel 3 -> UE1_Bs    channel 4 -> UE2_Bs   
    chOut=zeros(chTaps(1),length(waveform(:,1)));
    chBS_T0=zeros(chTaps(1),length(waveform(:,1)));
    TsVect=0:Ts:Ts*(length(waveform(:,1))-1);

    for channel=1:1:4
        chBS_T=chBS_T0; %reset
        for antenna=1:1:chTaps(1)
            for path=1:1:chTaps(3)

                inX=TsVect-chan(channel).delay(antenna,1,path,1);
                inY=interp1(TsVect,waveform(:,channel),inX,'pship');
                chBS_T(antenna,:)=inY*chan(channel).coeff(antenna,1,path,1)+chBS_T(antenna,:);%channel bs-terminal (only MP, no interf)

            end
        end
        chOut=chOut+chBS_T;%channel MP+interf (sum all obtained channels in one)
    end


    %% no BF

    %add noise

    chOut = awgn(chOut, snr, 'measured');

    chOut=transpose(chOut);

    title='Constellation without beamforming on UE1';
    bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,graph,title);
    [numbError,ratio]=biterr(in_t1,bits);
    fprintf('\nNumber of errors no beam 1: %d',numbError);
    fprintf('\nBER no beam 1: %f',ratio);
    ber1_nobf = [ber1_nobf, ratio]

    title='Constellation without beamforming on UE2';
    bits=OFDMDemod(ofdmMod,chOut(:,end),modOrder,graph,title);
    [numbError,ratio]=biterr(in_t2,bits);
    fprintf('\nNumber of errors no beam 2: %d',numbError);
    fprintf('\nBER no beam 2: %f',ratio);
    ber2_nobf = [ber2_nobf, ratio]

    %% algoritmo LMS 

    optimalWeight1 = LMSalgorithm(chOut,waveform(:,3),numArrayElements,Pars.lambda);  
    optimalWeight2 = LMSalgorithm(chOut,waveform(:,4),numArrayElements,Pars.lambda);   

    for i=1:length(waveform(:,3))
        y1(i)=chOut(i,:)*(optimalWeight1);
        y2(i)=chOut(i,:)*(optimalWeight2);
    end

    title='Constellation with LMS on UE1';
    bits=OFDMDemod(ofdmMod,transpose(y1),modOrder,graph,title);
    [numbError,ratio1]=biterr(in_t1,bits);
    fprintf('\nNumber of errors LMS1: %d',numbError);
    fprintf('\nBER LMS1: %f',ratio1);
    ber1_lms = [ber1_lms, ratio1]

    title='Constellation with LMS on UE2';
    bits=OFDMDemod(ofdmMod,transpose(y2),modOrder,graph,title);
    [numbError,ratio2]=biterr(in_t2,bits);
    fprintf('\nNumber of errors LMS2: %d',numbError);
    fprintf('\nBER LMS2: %f\n',ratio2);
    ber2_lms = [ber2_lms, ratio2]

    snr_vec = [snr_vec, snr]

    if ratio1 == 0 || ratio2 == 0
        break
    end
end

figure
ylim([0 1])


semilogy(snr_vec,ber1_lms,'-*',snr_vec,ber1_nobf,'-o',snr_vec,ber2_lms,'-*',snr_vec,ber2_nobf,'-o')
xlabel('SNR (dB)')
ylabel('BER')
legend('BER of UE1 with LMS BF','BER of UE1 without LMS BF','BER of UE2 with LMS BF','BER of UE2 without LMS BF')
%pause(1000)
%delete(pos);

%end







