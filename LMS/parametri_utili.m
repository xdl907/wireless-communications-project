function [BERr, BERm,SNRdifferences] = parametri_utili(referenceSignalWithoutNoise,referenceSignal, measuredSignal)
%PARAMETRI_UTILI calcolo BER e SNR
%input:
%   - referenceSignalWithoutNoise: signal before application of awgn
%   function
%   - referenceSignal: signal with awgn
%   - measuredSignal: signal with applied interferences and noise
%output
%   - BERr: BER for the reference signal
%   - BERm: BER for the measured signal
%   - SNRdifferences: effect of interference on SNR
    
    BERr = biterr(referenceSignal, referenceSignalWithoutNoise);
    BERm = biterr(measuredSignal, referenceSignalWithoutNoise);
    SNRWithoutInterference = snr(referenceSignalWithoutNoise, referenceSignal);
    SNRWithInterference = snr(referenceSignalWithoutNoise, referenceSignal);
    SNRdifferences = SNRWithoutInterference - SNRWithInterference;
end

