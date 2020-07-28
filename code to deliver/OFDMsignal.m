function [ofdmMod,waveform,in] = OFDMsignal(FFTLength, NumSymbols,modOrder)
% INPUT: FFTLength: lenght of the FFT 
%          NumSymbols: number of symbols
%          modOrder: order of QAM modulation
%          
% OUTPUT   ofdmMod: parameters of OFDM modulation
%          waveform: OFDM waveform generated
%          in: bits encoded in the OFDM signal 


% OFDM configuration:
ofdmMod = comm.OFDMModulator('FFTLength', FFTLength, ...
    'NumGuardBandCarriers', [1;1], ... 
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [floor(FFTLength/4)], ...
    'Windowing', false, ...
    'NumSymbols', NumSymbols, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [2;(FFTLength-1)]);

% input bit source:
in = randi([0 1], log2(modOrder)*NumSymbols*(FFTLength-4), 1);% FFTLength-4 -> I have to remove 2 symb for the guard bands and 2 for the pilots

dataInput = qammod(in, modOrder, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
ofdmInfo = info(ofdmMod);
ofdmSize = ofdmInfo.DataInputSize;
dataInput = reshape(dataInput, ofdmSize);

% waveform generation:
pilotInput = ones(2, NumSymbols, 1);
waveform = ofdmMod(dataInput, pilotInput);



end

