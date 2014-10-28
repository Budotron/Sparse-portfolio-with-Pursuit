function [gvx]=globalMinimumVariancePortfolio (Sigma, onesvec)
A=[2*Sigma onesvec;
    onesvec' 0];
zerosvec=zeros(length(onesvec), 1);
focs=[zerosvec;1];
gvx=A\focs;
gvx(length(gvx))=[];