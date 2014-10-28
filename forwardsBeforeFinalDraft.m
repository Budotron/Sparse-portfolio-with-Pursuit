S=getBestIndex(muvec)
    returnsTank=[min(muvec) min(FSMuvec)];
    allInds=[];
    S=getBestIndex(muvec)
    returnsTank=-inf;
    indexTank=[];
    expansionCount=1;
    for desiredSparsity=5 :5 :10
    expansionCount=expansionCount+1
        while nnz(indexTank)<desiredSparsity
        S=listInds(S, prices)
        [~, scols]=size(S)
        
        while scols<expansionCount
            S=listInds(S, prices)
        end
        keyboard
            for scan=1:length(S)
                scan
                tempS=S(scan, :);
                Q2=expandSigma(tempS, Sigma);
                smallOnesVec=ones(length(tempS), 1);
                gvx=globalMinimumVariancePortfolio(Q2, smallOnesVec);
                trialReturn=FSMuvec(tempS)'*gvx
                test1=(trialReturn>returnsTank);
                test2=trialReturn<1;
                test3 =-1<trialReturn;               
                if (test1==1 && test2==1 && test3==1)
                    returnsTank=trialReturn;
                    indexTank=unique([indexTank tempS]);
                end
                if scan==length(S)
                    break
                end
                
            end
            indexTank
            S=indexTank
        end
    end
    
    
    %%%latest
     for scan=1:length(S)
                    tempS=S(scan, :);
                    Q2=expandSigma(tempS, Sigma);
                    smallOnesVec=ones(length(tempS), 1);
                    gvx=globalMinimumVariancePortfolio(Q2, smallOnesVec);
                    trialReturn=muvec(tempS)'*gvx;
                    if trialReturn>returnsTank
                        returnsTank=trialReturn
                        indexTank=unique([indexTank tempS])
                    end
                    pause
                end
                S=indexTank
                
%%before latest%%%

%             S=getBestIndex(FSMuvec);
%             returnsTank=min(FSMuvec);
%             indexTank=[]
%             while nnz(S)<desiredSparsity
%                 allInds=listInds(S, prices);
%                 scan=0;
%                 for testInd=1:length(allInds)
%                     scan=scan+1;
%                     testInd=allInds(scan, :);
%                     M2=forwardSelectionM(testInd, :);
%                     Q2=expandSigma(testInd, forwardSelectionSigma);
%                     smallTestPortfolio=ComputePortfolio(M2, Q2, beta);
%                     fullPortfolio=zeros(prices, 1);
%                     fullPortfolio(testInd)=smallTestPortfolio;
%                     trialReturn=returnOnInvestment(...
%                         forwardSelectionPurchaseDayPrices,purchaseDayPrices,...
%                         fullPortfolio);
%                     test1=(trialReturn>returnsTank);
%                     test2=trialReturn<1;
%                     test3 =-1<trialReturn;
%                     if (test1==1 && test2==1 && test3==1)
%                         returnsTank=trialReturn;
%                         indexTank=unique([indexTank testInd])
%                     end
%                     
%                 end
%                 if nnz(indexTank<desiredSparsity)
%                     S=indexTank
%                 else
%                     break
%                 end
%             end