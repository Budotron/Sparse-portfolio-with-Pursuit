function [x]=useFormula (bigSig, onesvec, muvec, i)
A=[2*bigSig muvec' onesvec; muvec 0 0; onesvec' 0 0];
zerosvec=zeros(length(onesvec), 1);
b=[zerosvec; i; 1];
x=A\b;
x(length(x)-1:length(x))=[];


