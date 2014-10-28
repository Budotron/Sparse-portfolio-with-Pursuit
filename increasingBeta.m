%We want to determine if, in comparison to 1. the traditional formulation,
%2. the least squares regularization formulation, 3. the 1/N portfolio, if
%asset selection by pursuit algorithms 1. results in better perfoming
%portfolios (out of sample), 2. always select the same securities under
%similar economic conditions (stability), 3. select the same securities as
%those selected by L1 regularization/exact solution
%Question 1
%How do the returns compare to a 1/N strategy with time horizons one month,
%3 months and 1 year? 
%How do the returns compare to lasso (1 month, etc)?
%Pursuit algorithms induce sparsity. How do the returns of portfolios
%chosen by pursuit algorithms compare to other sparse selection algorithms?
%Question 2
%Are the same assets selected from year to year? What does this say about
%market conditions? 
%Question 3
%How much overlap is there between the assets selected by L1
%regualrization and each of the pursuit algorithms? What does this (lack
%of) overlap say about the algorithms?

clc 

warning('off', 'all');
warning;
load names; load datesSets; load datasets;

% oneOnNPortfolio=cell(1, length(datasets)-1);
% ExpectedOneOnNReturns=zeros(1, length(datasets)-1);
% realOneOnNReturns_month=zeros(1, length(datasets)-1);
% realOneOnNReturns_3months=zeros(1, length(datasets)-1);
% realOneOnNReturns_year=zeros(1, length(datasets)-1);
% oneOnNVar=zeros(1, length(datasets)-1);

% MarkowitzPortfolio=cell(1, length(datasets)-1);
% realMarkowitzReturns_month=zeros(1, length(datasets)-1);
% realMarkowitzReturns_3months=zeros(1, length(datasets)-1);
% realMarkowitzReturns_year=zeros(1, length(datasets)-1);
% MarkowitzVar=zeros(1, length(datasets)-1);
% 
% backwardsPort_all=cell(length(datasets)-1, length(datasets));
% BackwardsVar=zeros(length(datasets)-1, length(datasets));   
% realBackwardsReturns_month=zeros(length(datasets)-1, length(datasets));
% realBackwardsReturns_3months=zeros(length(datasets)-1, length(datasets));
% realBackwardsReturns_year=zeros(length(datasets)-1, length(datasets));
% backwardsSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% FowardsPort_all=cell(length(datasets)-1, length(datasets));
% ForwardsVar=zeros(length(datasets)-1, length(datasets));   
% realForwardsReturns_month=zeros(length(datasets)-1, length(datasets));
% realForwardsReturns_3months=zeros(length(datasets)-1, length(datasets));
% realForwardsReturns_year=zeros(length(datasets)-1, length(datasets));
% ForwardsSelectedAssets=cell(length(datasets)-1, length(datasets));

% OMPX_all=cell(length(datasets)-1, length(datasets));
% OMPVar=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_year=zeros(length(datasets)-1, length(datasets));
% OMPSelectedAssets=cell(length(datasets)-1, length(datasets));

% MPX_all=cell(length(datasets)-1, length(datasets));
% MPVar=zeros(length(datasets)-1, length(datasets));
% realMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realMPReturns_year=zeros(length(datasets)-1, length(datasets));
% MPSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% LSOMPX_all=cell(length(datasets)-1, length(datasets));
% LSOMPVar=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_year=zeros(length(datasets)-1, length(datasets));
% LSOMPSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% thrX_all=cell(length(datasets)-1, length(datasets));
% thrVar=zeros(length(datasets)-1, length(datasets));
% realthrReturns_month=zeros(length(datasets)-1, length(datasets));
% realthrReturns_3months=zeros(length(datasets)-1, length(datasets));
% realthrReturns_year=zeros(length(datasets)-1, length(datasets));
% thrSelectedAssets=cell(length(datasets)-1, length(datasets));

threshold=1e-8;



    thisPeriodPrices=datasets{1};
    thisPeriod=datesSets{1};
%     nextPeriodPrices=datasets{dataSetNumber+1};
%     nextPeriod=datesSets{dataSetNumber+1};
    [numPrices, numAssets]=size(thisPeriodPrices);
    investmentPrice=thisPeriodPrices(numPrices, :); %price today
    investmentDay=thisPeriod(numPrices, :);
%     rewardPrice_year=nextPeriodPrices(numPrices, :);
%     rewardDay_year=nextPeriod(numPrices, :);
%     rewardPrice_1month=nextPeriodPrices(30, :);
%     rewardDay_1month=nextPeriod(30, :);
%     rewardPrice_3months=nextPeriodPrices(90, :);
%     rewardDay_3months=nextPeriod(90, :);
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
    betacount=0
    for beta=0.01:0.01:0.04 %<?these are the expected returns. We want to see how
                            %the variance changes in response to increasing
                            %them
        betacount=betacount+1
        
        oneOnNPortfolio{betacount}=repmat(1/numAssets, [numAssets,1]) 
        oneOnNVar(betacount)=...
         calculatePortfolioVarience(oneOnNPortfolio{betacount}, Sigma)
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% calculate unconstrained Markowitz portfolios %%%%%%%%%%
        MarkowitzPortfolio{betacount}=ComputePortfolio(M, Sigma, beta); 

        %Do NOT need to calculate ANY OTHER expected returns! By design, any 
        %portfolio calculated from this point out is supposed to have an
        %expected return of beta!!!!!

        %%%%%%%%%%%%%%%%%% Markowitz portfolio variance %%%%%%%%%%%%%%%%%%%%%%%
        MarkowitzVar(betacount)=...
            calculatePortfolioVarience(MarkowitzPortfolio{betacount}, Sigma)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%% Backwards deletion strategy %%%%%%%%%%%%%%%%%%%%%%
        sparsity_number=4;
        sparseVec=[100, 50, 10]
        for desiredSparsity=sparseVec
            sparsity_number=sparsity_number-1;
            if sparsity_number==3
                nonSupport=[];
            end
            backwardsPort=abs(MarkowitzPortfolio{betacount});          
            while nnz(backwardsPort)>desiredSparsity
                forDeletion=deleteThisIndex(backwardsPort, nonSupport);
                nonSupport=[nonSupport, forDeletion]; %add the index to the non-support
                M1=restrictM(nonSupport, M); %M1 is M with rows in the set ns 
                                   %deleted  
                Q1=restrictSigma(Sigma, nonSupport); %Q1 is sigma with rows and 
                                           %cols ns deleted
                y=ComputePortfolio(M1, Q1, beta);%y is the portfolio 
                                               %obtained from ignoring
                                               %the contribution of 
                                               %assets in ns
                backwardsPort=reInsert_ns(y, nonSupport);%insert 0s into the non-support 
                                                %elements of y  
            end
            %%%%%%%%%%%%%%%%%% backwards portfolio variance %%%%%%%%%%%%%%%%%%%%
            BackwardsVar(betacount, 4-sparsity_number)=...
                        calculatePortfolioVarience(backwardsPort, Sigma);     
         
            end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        sparsity_number=0;
        sparsityVec=[10, 50, 100]
        for desiredSparsity=sparsityVec
            sparsity_number=sparsity_number+1
%             
%                     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     %%%%%%%%%%%%%%%%%%%% Forwards selection strategy %%%%%%%%%%%%%%%%%%
%                     if desiredSparsity==10
%                         S=getBestIndex(muvec);
%                     end
%                     while nnz(S)<desiredSparsity
%                         nnz(S);
%                         S=listInds(S, numAssets); %want to update this until 
%                                               %numrows(allInds)=desiredSparsity
%                         gmvVec=getGMVreturns(Sigma, S, muvec);
%                         [~, bestInd]=max(gmvVec);
%                         S=S(bestInd, :)         
%                         M2=M(S, :);
%                         Q2=expandSigma(S, Sigma);
%                         y=ComputePortfolio(M2, Q2, beta);
%                         ForwardsPort=buildForwardsX(y,S, numAssets);
%                     end
%                     %%%%%%%%%%%%%%%%%% forwards portfolio variance %%%%%%%%%%%%%%%%%%%%
%                     ForwardsVar(betacount, sparsity_number)=...
%                     calculatePortfolioVarience(ForwardsPort, Sigma);    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%% Pursuit selection setup%%%%%%%%%%%%%%%%%%%
        A=retmat;
        W=sqrt(diag(A'*A));
        for k=1:1:numAssets, 
            A(:,k)=A(:,k)/W(k); 
        end;
        b=beta*ones(numPrices-1, 1);
        OMPX=zeros(numAssets, 1);
        r=b;    
        %if desiredSparsity==10
            SS=[];%mark-about 3
        %end
        while nnz(OMPX)<desiredSparsity
            if  r'*r>threshold,
                Z=abs(A'*r);
                posZ=find(Z==max(Z));
                SS=sort([SS,posZ(1)]);
                r=b-A(:,SS)*pinv(A(:,SS))*b;    
            end;
            OMPX(SS)=pinv(A(:,SS))*b;
        end
        %%%%%%%%%%%%%%%%%%%%% Implement Part (2) %%%%%%%%%%%%%%%%%%%%%%
        M3=M(SS, :);
        Q3=expandSigma(SS, Sigma);
        y=ComputePortfolio(M3, Q3, beta);
        OMPX(SS)=y;
        selectedAssets=find(OMPX);
        OMPVar(betacount, sparsity_number)=...
                calculatePortfolioVarience(OMPX, Sigma) 
                %%%%%%%%%%%%%%%%%%%%%%%% MP selection strategy %%%%%%%%%%%%%%%%%%%%%
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
        MPVar(betacount, sparsity_number)=...
                calculatePortfolioVarience(MPX, Sigma);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
        %%%%%%%%%%%%%%%%%%%%%% LSOMP selection strategy %%%%%%%%%%%%%%%%%%%
        
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
        LSOMPVar(betacount, sparsity_number)=...
            calculatePortfolioVarience(LSOMPX, Sigma);
        selectedAssets=find(LSOMPX);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                       

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
           
            thrVar(betacount, sparsity_number)=...
            calculatePortfolioVarience(thrX, Sigma);
                     
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
% 
    end
