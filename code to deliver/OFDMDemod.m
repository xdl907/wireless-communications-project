function [dataOutput] = OFDMDemod(ofdmMod,waveform,modOrder,scatterPlot,title)
% INPUT: ofdmMod: parameters of the OFDM modulation
%          waveform: samples of the signal to be demodulated 
%          modOrder: order of QAM modulation
%          scatterPlot: boolean, if true the function shows the
%          constellation plot
%          title: title of the constellation plot
%          

% OUTPUT   dataOutput: bits contained in the waveform 


ofdmDemod=comm.OFDMDemodulator(ofdmMod);

%extract symbols from OFDM signal
symbols=ofdmDemod(waveform);

x=real(symbols);
x=reshape(x,[],1);
y=imag(symbols);
y=reshape(y,[],1);

if scatterPlot==true
    figure('Name',title);
    scatter(x,y,100,'r','LineWidth',0.01)
    
end

dataOutput = qamdemod(symbols, modOrder, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true);
dataOutput = reshape(dataOutput, [],1);

end

