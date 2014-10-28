clc, 
clear all
K=[50, 100, 150];
kcount=0;
for k=K
    kcount=kcount+1
    load datapoints
    [~,totalAssets]=size(datapoints);
    choose=randsample(1:totalAssets, k)
    load dataSets
    
    oneOnNPortfolio=cell(length(K), length(datasets)-1);
    ExpectedOneOnNReturns=zeros(length(K), length(datasets)-1);
    realOneOnNReturns_month=zeros(length(K), length(datasets)-1);
    realOneOnNReturns_3months=zeros(length(K), length(datasets)-1);
    realOneOnNReturns_year=zeros(length(K), length(datasets)-1);
    oneOnNVar=zeros(length(K), length(datasets)-1);

    MarkowitzPortfolio=cell(length(K), length(datasets)-1);
    realMarkowitzReturns_month=zeros(length(K), length(datasets)-1);
    realMarkowitzReturns_3months=zeros(length(K), length(datasets)-1);
    realMarkowitzReturns_year=zeros(length(K), length(datasets)-1);
    MarkowitzVar=zeros(length(K), length(datasets)-1);
    
    FowardsPort_all=cell(length(datasets)-1, length(datasets));
    FowardsPort_allKays=cell(length(K), length(datasets));
    ForwardsVar=zeros(length(datasets)-1, length(datasets));   
    realForwardsReturns_month=zeros(length(datasets)-1, length(datasets));
    realForwardsReturns_3months=zeros(length(datasets)-1, length(datasets));
    realForwardsReturns_year=zeros(length(datasets)-1, length(datasets));
    ForwardsSelectedAssets=cell(length(datasets)-1, length(datasets));
    
    threshold=1e-8;
    dataSetCount=0;
    for dataSetNumber=1:length(datasets)-1
        dataSetCount=dataSetCount+1
        thisPeriodPrices=datasets{dataSetNumber}(:, choose);
        nextPeriodPrices=datasets{dataSetNumber+1}(:, choose);
        [numPrices, numAssets]=size(thisPeriodPrices);
        investmentPrice=thisPeriodPrices(numPrices, :); %price today
        rewardPrice_year=nextPeriodPrices(numPrices, :);
        rewardPrice_1month=nextPeriodPrices(30, :);
        rewardPrice_3months=nextPeriodPrices(90, :); 
        retmat=makeReturnsMatrix(thisPeriodPrices);    
        muvec=mean(retmat)'; %vector of expected returns
        Sigma=cov(retmat);   %covariance matrix of return vectors
        onesvec=ones(length(muvec), 1);    
        gvx=globalMinimumVariancePortfolio(Sigma, onesvec); %portfolio from 
                                                            %ignoring
                                                            %target constraint
        betaMin=muvec'*gvx; %expected return of GMV portfolio. Any portfolio 
                            %with a higer expected return will be non-dominated
        M=[muvec onesvec]; %concatenation for use in formula     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%   set up 1/N portfolio    %%%%%%%%%%%%%%%%%%%%%%%
        oneOnNPortfolio{kcount,dataSetNumber}=repmat(1/numAssets, [numAssets,1]);  
        %%%%%%%%%%%%%%% expected returns on 1/N ports per period %%%%%%%%%%%%%
        ExpectedOneOnNReturns(kcount, dataSetNumber)=muvec'*oneOnNPortfolio{kcount,dataSetNumber} ;    

        %%%%%%%%%%%%%%%%%%%%%%% 1/N portfolio variance %%%%%%%%%%%%%%%%%%%%%%%%
        oneOnNVar(kcount, dataSetNumber)=...
             calculatePortfolioVarience(oneOnNPortfolio{kcount,dataSetNumber}, Sigma);

        %%%%%%%%%%%%%%%% calculate returns on 1/N ports per period %%%%%%%%%%%%
        realOneOnNReturns_month(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_1month,...
        oneOnNPortfolio{kcount, dataSetNumber});
        realOneOnNReturns_3months(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_3months,...
        oneOnNPortfolio{kcount, dataSetNumber});
        realOneOnNReturns_year(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_year,...
        oneOnNPortfolio{kcount, dataSetNumber});

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        beta=ExpectedOneOnNReturns(dataSetNumber)  ;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%% calculate unconstrained Markowitz portfolios %%%%%%%%%%
        MarkowitzPortfolio{kcount, dataSetNumber}=ComputePortfolio(M, Sigma, beta); 

        %Do NOT need to calculate ANY OTHER expected returns! By design, any 
        %portfolio calculated from this point out is supposed to have an
        %expected return of beta!!!!!

        %%%%%%%%%%%%%%%%% Markowitz portfolio variance %%%%%%%%%%%%%%%%%%%%%%%
        MarkowitzVar(kcount, dataSetNumber)=...
            calculatePortfolioVarience(MarkowitzPortfolio{kcount, dataSetNumber}, Sigma);
        %%%%%%%%%%%%%%% calculate returns on Markowitz portfolios %%%%%%%%%%%%
        realMarkowitzReturns_month(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_1month,...
        MarkowitzPortfolio{kcount, dataSetNumber});
        realMarkowitzReturns_3months(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_3months,...
        MarkowitzPortfolio{kcount, dataSetNumber});
        realMarkowitzReturns_year(kcount, dataSetNumber)=...
        returnOnInvestment(investmentPrice, rewardPrice_year,...
        MarkowitzPortfolio{kcount, dataSetNumber});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%% Backwards deletion strategy %%%%%%%%%%%%%%%%%%%%%%
%         sparsity_number=11;
%         for desiredSparsity=k/2:-10:k/10
%             sparsity_number=sparsity_number-1
%             if sparsity_number==10
%                 nonSupport=[];
%             end
%             backwardsPort=abs(MarkowitzPortfolio{dataSetNumber});          
%             while nnz(backwardsPort)>desiredSparsity
%                 forDeletion=deleteThisIndex(backwardsPort, nonSupport);
%                 nonSupport=[nonSupport, forDeletion]; %add the index to the non-support
%                 M1=restrictM(nonSupport, M); %M1 is M with rows in the set ns deleted  
%                 Q1=restrictSigma(Sigma, nonSupport); %Q1 is sigma with rows and  cols ns deleted
%                 y=ComputePortfolio(M1, Q1, beta);
%                 %y is the portfolio obtained from ignoring the contribution of 
%                                               %assets in ns
%                 backwardsPort=reInsert_ns(y, nonSupport);%insert 0s into the non-support 
%                                                % elements of y  
%             end
%             backwardsPort,
%             nnz(backwardsPort), keyboard
%             %backwardsPort_all{dataSetNumber, 11-sparsity_number}=backwardsPort;
%             selectedAssets=find(backwardsPort);
%             backwardsSelectedAssets{dataSetNumber, 11-sparsity_number}=...
%                 names(selectedAssets);
%             %%%%%%%%%%%%%%%%% backwards portfolio variance %%%%%%%%%%%%%%%%%%%%
%             BackwardsVar(dataSetNumber, 11-sparsity_number)=...
%                       calculatePortfolioVarience(backwardsPort, Sigma);     
%             %%%%%%%%%%%%%%%%% backwards portfolio returns %%%%%%%%%%%%%%%%%%%%
%             realBackwardsReturns_month(dataSetNumber, 11-sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_1month,...
%             backwardsPort);
%             realBackwardsReturns_3months(dataSetNumber, 11-sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_3months,...
%             backwardsPort);
%             realBackwardsReturns_year(dataSetNumber, 11-sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_year,...
%             backwardsPort);          
%         end
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        sparsity_number=0;
        for desiredSparsity=k/10:5:k/5 
            sparsity_number=sparsity_number+1
        
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% BNB solution by YALMIP %%%%%%%%%%%%%%%%%%%%%%
            BNBX=branchAndBoundSolve(Sigma, beta, muvec', desiredSparsity);
            capture(dataSetNumber, sparsity_number)=nnz(BNBX);
            assetsChosenBNBX{dataSetNumber, sparsity_number}=names(find(BNBX))
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
            %%%%%%%%%%%%%%%%%% Forwards selection strategy %%%%%%%%%%%%%%%%%%
            if desiredSparsity==k/10
                S=getBestIndex(muvec);
            end
            while nnz(S)<desiredSparsity
                %nnz(S);
                S=listInds(S, numAssets); %want to update this until 
                                      %numrows(allInds)=desiredSparsity
                gmvVec=getGMVreturns(Sigma, S, muvec);
                [~, bestInd]=max(gmvVec);
                S=S(bestInd, :);         
                M2=M(S, :);
                Q2=expandSigma(S, Sigma);
                y=ComputePortfolio(M2, Q2, beta);
                ForwardsPort=buildForwardsX(y,S, numAssets);
            end
            ForwardsPort_all{dataSetNumber, sparsity_number}=ForwardsPort;
            selectedAssets=find(ForwardsPort);
            ForwardsSelectedAssets{dataSetNumber, sparsity_number}=...
                names(selectedAssets);
            %%%%%%%%%%%%%%%% forwards portfolio variance %%%%%%%%%%%%%%%%%%%%
            ForwardsVar(dataSetNumber, sparsity_number)=...
                        calculatePortfolioVarience(ForwardsPort, Sigma);     
            %%%%%%%%%%%%%%%% forwards portfolio returns %%%%%%%%%%%%%%%%%%%%
            realForwardsReturns_month(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_1month,...
            ForwardsPort);
            realForwardsReturns_3months(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_3months,...
            ForwardsPort);
            realForwardsReturns_year(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_year,...
            ForwardsPort);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%% Pursuit selection setup%%%%%%%%%%%%%%%%%%%
            A=retmat;
            W=sqrt(diag(A'*A));
            for k=1:1:numAssets, 
                A(:,k)=A(:,k)/W(k); 
            end;
            b=beta*ones(numPrices-1, 1);

%             %%%%%%%%%%%%%%%%%%% all selection strategy %%%%%%%%%%%%%%%%%%%
%             idea:            
%             for each dataset
%                for each target return
%                    for each desired sparsity
%                        use PA to select a K-sparse approximant to the vector 
%                        b(=beta) with respect to the dictionary formed by the 
%                        matrix of all returns(=retmat)
%                    end
%                    use formula to reclaibrate the K-sparse approximant into
%                    a portfolio with target return beta
%                end
%                calculate the returns over the time horizons
%             end
%             %%%%%%%%%%%%%%%%%%%% Implement Part (1) %%%%%%%%%%%%%%%%%%%%%%
            OMPX=zeros(numAssets, 1);
            r=b;    
            if desiredSparsity==k/10
                SS=[];%mark-about 3
            end
            while nnz(OMPX)<desiredSparsity
                if  r'*r>threshold,
                    Z=abs(A'*r);
                    posZ=find(Z==max(Z));
                    SS=sort([SS,posZ(1)]);
                    r=b-A(:,SS)*pinv(A(:,SS))*b;    
                end;
                OMPX(SS)=pinv(A(:,SS))*b;
            end
            %%%%%%%%%%%%%%%%%%%% Implement Part (2) %%%%%%%%%%%%%%%%%%%%%%
            M3=M(SS, :);
            Q3=expandSigma(SS, Sigma);
            y=ComputePortfolio(M3, Q3, beta);
            OMPX(SS)=y;
            OMPX_all{dataSetNumber, sparsity_number}=OMPX;
            selectedAssets=find(OMPX);
            OMPSelectedAssets{dataSetNumber, sparsity_number}=...
                names(selectedAssets);
            %%%%%%%%%%%%%%%%%%%% Implement Part (3) %%%%%%%%%%%%%%%%%%%%%%
            OMPVar(dataSetNumber, sparsity_number)=...
                    calculatePortfolioVarience(OMPX, Sigma);     
            %%%%%%%%%%%%%%%%% OMP portfolio returns %%%%%%%%%%%%%%%%%%%%
            realOMPReturns_month(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_1month,...
            OMPX);
            realOMPReturns_3months(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_3months,...
            OMPX);
            realOMPReturns_year(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_year,...
            OMPX);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            %%%%%%%%%%%%%%%%%%%%%%% MP selection strategy %%%%%%%%%%%%%%%%%%%%%
            MPX=zeros(numAssets,1);
            r=b;
            while nnz(MPX)<desiredSparsity
                if r'*r>threshold,
                    Z=abs(A'*r);
                    posZ=find(Z==max(Z),1);
                    MPX(posZ)=MPX(posZ)+A(:,posZ)'*r;
                    r=r-A(:,posZ)*A(:,posZ)'*r;
                end;
                SS=find(abs(MPX)>1e-8)';
            end
            M4=M(SS, :);
            Q4=expandSigma(SS, Sigma);
            y=ComputePortfolio(M4, Q4, beta);
            MPX(SS)=y; 
            MPX_all{dataSetNumber, sparsity_number}=MPX;
            selectedAssets=find(MPX);
            MPVar(dataSetNumber, sparsity_number)=...
                    calculatePortfolioVarience(MPX, Sigma);
            MPSelectedAssets{dataSetNumber, sparsity_number}=...
                names(selectedAssets);     
            %%%%%%%%%%%%%%%%% MP portfolio returns %%%%%%%%%%%%%%%%%%%%
            realMPReturns_month(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_1month,...
            MPX);
            realMPReturns_3months(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_3months,...
            MPX);
            realMPReturns_year(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_year,...
            MPX);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%% LSOMP selection strategy %%%%%%%%%%%%%%%%%%%

            r=b;
            if sparsity_number==1 
                SS=[];
                nnz(SS)
            end
            while nnz(SS)<desiredSparsity
                if r'*r>threshold,
                    Z=zeros(numAssets,1);
                    for jj=1:1:numAssets
                        SStemp=[SS,jj];
                        rtemp=b-A(:,SStemp)*pinv(A(:,SStemp))*b;
                        Z(jj)=rtemp'*rtemp;
                    end;
                    posZ=find(Z==min(Z),1);
                    SS=sort([SS,posZ(1)]);
                    r=b-A(:,SS)*pinv(A(:,SS))*b;    
                end;
                LSOMPX=zeros(numAssets, 1);
                LSOMPX(SS)=pinv(A(:,SS))*b;
            end
            nnz(SS)
            M5=M(SS, :);
            Q5=expandSigma(SS, Sigma);
            y=ComputePortfolio(M5, Q5, beta);
            LSOMPX(SS)=y; 
            LSOMPX_all{dataSetNumber, sparsity_number}=LSOMPX;
            LSOMPVar(dataSetNumber, sparsity_number)=...
                calculatePortfolioVarience(LSOMPX, Sigma);
            selectedAssets=find(LSOMPX);
            LSOMPSelectedAssets{dataSetNumber, sparsity_number}=...
                names(selectedAssets)   ;
            %%%%%%%%%%%%%%%%%%% LSOMP portfolio returns %%%%%%%%%%%%%%%%%%%%%%
            realLSOMPReturns_month(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_1month,...
            LSOMPX);
            realLSOMPReturns_3months(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_3months,...
            LSOMPX);
            realLSOMPReturns_year(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_year,...
            LSOMPX);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%% tresholding selection strategy %%%%%%%%%%%%%%
            thrX=zeros(numAssets, 1); 
            Z=A'*b; 
            [Za,posZ]=sort(abs(Z),'descend');
            r=b;
            SS=[];
            while nnz(thrX)<desiredSparsity
                if r'*r>threshold, 
                    SS=[SS,posZ(length(SS)+1)];
                    thrX(SS)=pinv(A(:,SS))*b;
                    r=b-A(:,SS)*thrX(SS);
                end;                
            end
            M6=M(SS, :);
            Q6=expandSigma(SS, Sigma);
            y=ComputePortfolio(M6, Q6, beta);
            thrX(SS)=y;
            thrX_all{dataSetNumber, sparsity_number}=thrX;
            thrVar(dataSetNumber, sparsity_number)=...
            calculatePortfolioVarience(thrX, Sigma);
            selectedAssets=find(thrX);
            thrSelectedAssets{dataSetNumber, sparsity_number}=...
                names(selectedAssets)   ;
            %%%%%%%%%%%%%%%%%%%% thr portfolio returns %%%%%%%%%%%%%%%%%%%%%%
            realthrReturns_month(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_1month,...
            thrX);
            realthrReturns_3months(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_3months,...
            thrX);
            realthrReturns_year(dataSetNumber, sparsity_number)=...
            returnOnInvestment(investmentPrice, rewardPrice_year,...
            thrX);           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
         keyboard
    end
end
    