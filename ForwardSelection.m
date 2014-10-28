% idea: 
% for each dataset
%     for each target return
%         for each desired sparsity
%             find the best perfoming global minimum variance portfolio \leq
%             desired sparsity
%         end
%         use formula to recalibrate weights so that the return on the portfolio 
%         matches target return
%     end
%     calculate the return over the time horizons
% end
            
clc, clear all
tempData=cell(3, 1);
for i=1:length(tempData)
    tempData{i}=Generatedata(10, 18);
end
datasets=tempData;
for dataSetNumber=1:length(tempData)
    dailyPrices=datasets{dataSetNumber};
    [numDays,prices]=size(dailyPrices);
    lastTenth=1/10*numDays;
    purchaseDayPrices=dailyPrices(numDays-lastTenth, :);
    retmat=makeReturnsMatrix(dailyPrices(1:numDays-lastTenth, :));
    Sigma=cov(retmat);
    muvec=mean(retmat)';
    onesvec=ones(length(muvec), 1);
    M=[muvec onesvec];
    gvx=globalMinimumVariancePortfolio(Sigma, onesvec);
    betaMin=muvec'*gvx;
    betaMax=0.1;
    betarange=sort([betaMin, betaMax], 'ascend');
    targetReturnNumber=0;
    for beta=min(betarange):(max(betarange)-min(betarange))/2:max(betarange)
        targetReturnNumber=targetReturnNumber+1;
        sparsity_set=0;
        expansionCount=1;
        for desiredSparsity=5:5:15
            sparsity_set=sparsity_set+1;
            %generate candidate indices
            S=getBestIndex(muvec);
            indexTank=[];
            expansionCount=expansionCount+1
            while nnz(S)<desiredSparsity
                S=listInds(S, prices) %want to update this until 
                                      %numrows(allInds)=desiredSparsity
                
                                      
                gmvVec=getGMVreturns(Sigma, S, muvec)
                pause
                [~, bestInd]=max(gmvVec);
                S=S(bestInd, :);
            end
            solutionSparsity=nnz(unique(S));
            [desiredSparsity solutionSparsity]
        end
    end
end
