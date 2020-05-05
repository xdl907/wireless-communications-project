function [dataOutput] = OFDMDemod(ofdmMod,waveform,modOrder,scatterPlot,title)

ofdmDemod=comm.OFDMDemodulator(ofdmMod);

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

