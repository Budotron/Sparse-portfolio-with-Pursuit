clc 
%clear all
warning('off', 'all');
warning;
load names; load datesSets; load datasets;
assetsChosenBNBX=cell(9, 10)
capture=(zeros(10, 10));
for dataSetNumber=1:length(datasets)-1
    tic
%     thisPeriodPrices=datasets{dataSetNumber};
%     thisPeriod=datesSets{dataSetNumber};
%     [numPrices, numAssets]=size(thisPeriodPrices);
numAssets=100;
    thisPeriodPrices=randn(numAssets)
    retmat=makeReturnsMatrix(thisPeriodPrices) ;  
    muvec=mean(retmat)'; %vector of expected returns
    Sigma=cov(retmat);   %covariance matrix of return vectors
    onesvec=ones(length(muvec), 1); 
    gvx=globalMinimumVariancePortfolio(Sigma, onesvec); %portfolio from 
                                                        %ignoring
                                                        %target constraint
    betaMin=muvec'*gvx %expected return of GMV portfolio. Any portfolio 
                        %with a higer expected return will be non-dominated  
    M=[muvec onesvec]; %concatenation for use in formula 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   set up 1/N portfolio    %%%%%%%%%%%%%%%%%%%%%%%
    oneOnNPortfolio{dataSetNumber}=repmat(1/numAssets, [numAssets,1])  
    
    %%%%%%%%%%%%%%%% expected returns on 1/N ports per period %%%%%%%%%%%%%
    ExpectedOneOnNReturns(dataSetNumber)=muvec'*oneOnNPortfolio{dataSetNumber}
    beta=ExpectedOneOnNReturns(dataSetNumber)
    sparsity_number=0;
    for desiredSparsity=4:2:8
        sparsity_number=sparsity_number+1
        BNBX=branchAndBoundSolve(Sigma, beta, muvec', desiredSparsity);
        capture(dataSetNumber, sparsity_number)=nnz(BNBX);
        assetsChosenBNBX{dataSetNumber, sparsity_number}=names(find(BNBX))
        
    end
    toc
end
save assetsChosenBNBX
save capture
