function [waveform,ofdmMod, Fs] = OFDMsignal()
% Generated by MATLAB(R) 9.7 (R2019b) and Communications Toolbox 7.2 (R2019b).
% Generated on: 08-Mar-2020 16:25:04

%% Generate OFDM Waveform
% OFDM configuration:
ofdmMod = comm.OFDMModulator('FFTLength', 64, ...
    'NumGuardBandCarriers', [6;5], ...
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [16], ...
    'Windowing', false, ...
    'NumSymbols', 100, ...
    'NumTransmitAntennas', 4, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [12;26;40;54]);

M = 4; 	 % Modulation order
% input bit source:
in = randi([0 1], 9800, 4);

dataInput = qammod(in, M, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
ofdmInfo = info(ofdmMod);
ofdmSize = ofdmInfo.DataInputSize;
dataInput = reshape(dataInput, ofdmSize);

% waveform generation:
pilotInput = ones(4, 100, 4);
waveform = ofdmMod(dataInput, pilotInput);

Fs = 6.4e+07; 								 % sample rate of waveform



end

