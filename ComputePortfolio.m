function[x]=ComputePortfolio(M, Sigma, beta)
mutilde=[beta;1];
x=inv(Sigma)*M*inv(M'*inv(Sigma)*M)*mutilde;