function [forwardPortfolio]=buildForwardsX(tempPort, testSupport, ...
                prices)
            
x=zeros(prices, 1);
x(testSupport)=tempPort;
forwardPortfolio=x;
