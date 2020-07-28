function [w] = LMSalgorithm(chOut,r,numArrayElements,~)


w=zeros(numArrayElements*numArrayElements,1);


 
for i=1:length(r)
    

    y=chOut(i,:)*w;
    e=r(i)-y;
mu=2/(1000*max(eig(chOut(i,:)*chOut(i,:)')));

    w=w+mu*chOut(i,:)'*e;
    

end


end

