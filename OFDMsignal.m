function [ofdmMod,waveform,in] = OFDMsignal(FFTLength, NumSymbols,modOrder)
% Generated by MATLAB(R) 9.7 (R2019b) and Communications Toolbox 7.2 (R2019b).
% Generated on: 08-Mar-2020 16:25:04


% OFDM configuration:
ofdmMod = comm.OFDMModulator('FFTLength', FFTLength, ...
    'NumGuardBandCarriers', [1;1], ... %1 guard band all'inizio e una alla fine (FFTLength)
    'InsertDCNull', false, ...
    'CyclicPrefixLength', [floor(FFTLength/2)], ...
    'Windowing', false, ...
    'NumSymbols', NumSymbols, ...
    'NumTransmitAntennas', 1, ...
    'PilotInputPort', true, ...
    'PilotCarrierIndices', [2;(FFTLength-1)]);%un segnale di pilot appena dopo la prima guard band (quindi in pos 2) e uno nel penultimo simbolo (FFTLength-1)

% input bit source:
in = randi([0 1], log2(modOrder)*NumSymbols*(FFTLength-4), 1);% FFTLength-4 >> devo togliere i 2 simboli per le guard bands e i 2 per i pilot

dataInput = qammod(in, modOrder, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
ofdmInfo = info(ofdmMod);
ofdmSize = ofdmInfo.DataInputSize;
dataInput = reshape(dataInput, ofdmSize);

% waveform generation:
pilotInput = ones(2, NumSymbols, 1);
waveform = ofdmMod(dataInput, pilotInput);



end
