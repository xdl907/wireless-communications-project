function [w] = LMS2(chOut,r, numArrayElements,lambda)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

w=zeros(numArrayElements*numArrayElements,1);

chOut=transpose(chOut);

%mu = (2/(lambda))/10000;
 
for i=1:length(r)
    
    y=w'*chOut(:,i);
    e=r(i)-y;
    mu=2/(1000*max(eig(chOut(:,i)*chOut(:,i)')));
   

    w=w+mu*(chOut(:,i))*conj(e);

end

end

