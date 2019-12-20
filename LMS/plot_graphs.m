function plot_graphs(array, carrierFrequency, elevation,weights, arrivingSignal, stimatedSignal)
%PLOT_GRAPHS make the plots for the input data for beamforming
%input:
%   - array: antenna array for the systen
%   - carrierFrequency
%   - elevation: of the terminal wrt BS
%   - weights: for the plot of beam
%   - arrivingSignal: with interferences and noise
%   - stimatedSignal: only with noises
%output: plotting of the graphs
tiledlayout(2,1);
nexttile
pattern(array,carrierFrequency,[-180:180],elevation,'PropagationSpeed',carrierFrequency,'Type',...
    'powerdb','CoordinateSystem','polar','Weights',weights)
title('beamgraph of BS-terminal')
% plot of the signal
nexttile
plot(real(arrivingSignal));
hold on;
plot(real(stimatedSignal));
title('comparisons of arrived signal and stimated signal')
end

