function [gmvVec]=getGMVreturns(Sigma, S, muvec)
  for scan=1:length(S)
      tempS=S(scan, :);
      Q2=expandSigma(tempS, Sigma);
      smallOnesVec=ones(length(tempS), 1);
      gmv=globalMinimumVariancePortfolio(Q2, smallOnesVec);
      gmvVec(scan)=muvec(tempS)'*gmv;
      hold=bopRepeats(S);
      gmvVec(hold)=-inf;
  end
