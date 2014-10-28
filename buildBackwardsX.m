function [backwardsPortfolio]=buildBackwardsX(tempPort, testSupport, ...
                MarkowitzPortfoliosMat, dataSet, targetReturnNumber)
            
x=zeros(length(MarkowitzPortfoliosMat{dataSet,targetReturnNumber}), 1);
x(testSupport)=tempPort;
backwardsPortfolio=x
