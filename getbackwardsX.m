function [feedforward]= getbackwardsX(x, ns, desiredSparsity)
while nnz(x)>desiredSparsity
                  %find the index of the smallest asset weight in the 
                  %portfolio
                  forDeletion=deleteThisIndex(x, ns);
                  ns=[ns, forDeletion]; %add the index to the non-support
                  M1=restrictM(ns, M); %M1 is M with rows in the set ns 
                                       %deleted  
                  Q1=restrictSigma(Sigma, ns); %Q1 is sigma with rows and 
                                               %cols ns deleted
                  y=ComputePortfolio(M1, Q1, beta);%y is the portfolio 
                                                   %obtained from ignoring
                                                   %the contribution of 
                                                   %assets in ns
                  backwardsX=reInsert_ns(y, ns);%insert 0s into the non-support elements of y
                  feedforward=[backwardsX', ns]
end                                       
