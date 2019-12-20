function [wNew] = LMSalgorithm(chOut,r, wOld)
%LMSalgorithm: compute weights set
%input:
%   -chOut: channel output, a matrix where each column is the rx signal at
%   the ith array element
%   -r: reference signal
%output
%   -wNew: the updated array of weights

Y = chOut * transpose(conj(chOut));
R = corr(Y); % autocorrelation matrix of the rx signal at each array element
z = corr(chOut, transpose(r)); %correlation between the rx signals and the reference one 

% MSE gradient
MSEgrad = (2.*R*wOld)-(2.*z);

% step size mu
mu = 1 / trace(R);

%update weights 
wNew=wOld-0.5*mu*MSEgrad;

end

