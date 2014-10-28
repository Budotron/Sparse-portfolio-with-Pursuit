function [betaMin]=getBetaMin(bigSig, onesvec, muvec)
A=[2*bigSig onesvec; onesvec' 0];
zerosvec=zeros(length(onesvec), 1);
b=[zerosvec; 1];
x=A\b;
x(length(x))=[];
betaMin=muvec*x;

