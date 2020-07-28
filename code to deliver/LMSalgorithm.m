function [w] = LMSalgorithm(chOut,r,numArrayElements)
% INPUT: chOut: matrix of all samples received during the transmission by
%              all array elements
%          r: samples of the signal transmitted 
%          numArrayElements: Number of elements that compose the URA
%                            (size)
%          

% OUTPUT   w: weights

w=zeros(numArrayElements*numArrayElements,1); %initialize vector of weights


 
for i=1:length(r)
    

    y=chOut(i,:)*w; %weighted sample
    e=r(i)-y; %compute the error
    
    mu=2/(1000*max(eig(chOut(i,:)'*chOut(i,:)))); %step size

    w=w+mu*chOut(i,:)'*e; %update weights 
    

end


end

