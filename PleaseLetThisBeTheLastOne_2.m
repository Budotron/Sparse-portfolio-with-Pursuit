%We want to determine if, in comparison to 1. the traditional formulation,
%2. the least squares regularization formulation, 3. the 1/N portfolio, if
%asset selection by pursuit algorithms 1. results in better perfoming
%portfolios (out of sample), 2. always select the same securities under
%similar economic conditions (stability), 3. select the same securities as
%those selected by L1 regularization. 
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
%clear all
warning('off', 'all');
warning;
load names; load datesSets; load datasets;
realOneOnNReturns_month=zeros(1, length(datasets)-1);
realOneOnNReturns_3months=zeros(1, length(datasets)-1);
realOneOnNReturns_year=zeros(1, length(datasets)-1);
oneOnNVar=zeros(1, length(datasets)-1);

% realMarkowitzReturns_month=zeros(1, length(datasets)-1);
% realMarkowitzReturns_3months=zeros(1, length(datasets)-1);
% realMarkowitzReturns_year=zeros(1, length(datasets)-1);
% MarkowitzVar=zeros(1, length(datasets)-1);
% 
% BackwardsVar=zeros(length(datasets)-1, length(datasets));   
% realBackwardsReturns_month=zeros(length(datasets)-1, length(datasets));
% realBackwardsReturns_3months=zeros(length(datasets)-1, length(datasets));
% realBackwardsReturns_year=zeros(length(datasets)-1, length(datasets));
% backwardsSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% ForwardsVar=zeros(length(datasets)-1, length(datasets));   
% realForwardsReturns_month=zeros(length(datasets)-1, length(datasets));
% realForwardsReturns_3months=zeros(length(datasets)-1, length(datasets));
% realForwardsReturns_year=zeros(length(datasets)-1, length(datasets));
% ForwardsSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% OMPVar=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realOMPReturns_year=zeros(length(datasets)-1, length(datasets));
% OMPSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% MPVar=zeros(length(datasets)-1, length(datasets));
% realMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realMPReturns_year=zeros(length(datasets)-1, length(datasets));
% MPSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% LSOMPVar=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_month=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_3months=zeros(length(datasets)-1, length(datasets));
% realLSOMPReturns_year=zeros(length(datasets)-1, length(datasets));
% LSOMPSelectedAssets=cell(length(datasets)-1, length(datasets));
% 
% thrVar=zeros(length(datasets)-1, length(datasets));
% realthrReturns_month=zeros(length(datasets)-1, length(datasets));
% realthrReturns_3months=zeros(length(datasets)-1, length(datasets));
% realthrReturns_year=zeros(length(datasets)-1, length(datasets));
% thrSelectedAssets=cell(length(datasets)-1, length(datasets));

ExactVar=zeros(length(datasets)-1, length(datasets));
ExactReturns_month=zeros(length(datasets)-1, length(datasets));
ExactReturns_3months=zeros(length(datasets)-1, length(datasets));
exactReturns_year=zeros(length(datasets)-1, length(datasets));

threshold=1e-8;


for dataSetNumber=1:length(datasets)-1
    thisPeriodPrices=datasets{dataSetNumber};
    nextPeriodPrices=datasets{dataSetNumber+1};
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   set up 1/N portfolio    %%%%%%%%%%%%%%%%%%%%%%%
    oneOnNPortfolio=repmat(1/numAssets, [numAssets,1]);  
    
    %%%%%%%%%%%%%%%% expected returns on 1/N ports per period %%%%%%%%%%%%%
    ExpectedOneOnNReturns=muvec'*oneOnNPortfolio;      
    
    %%%%%%%%%%%%%%%%%%%%%%% 1/N portfolio variance %%%%%%%%%%%%%%%%%%%%%%%%
    oneOnNVar(dataSetNumber)=...
                        calculatePortfolioVarience(oneOnNPortfolio, Sigma);
    
    %%%%%%%%%%%%%%%% calculate returns on 1/N ports per period %%%%%%%%%%%%
    realOneOnNReturns_month(dataSetNumber)=...
    returnOnInvestment(investmentPrice, rewardPrice_1month,...
    oneOnNPortfolio);
    realOneOnNReturns_3months(dataSetNumber)=...
    returnOnInvestment(investmentPrice, rewardPrice_3months,...
    oneOnNPortfolio);
    realOneOnNReturns_year(dataSetNumber)=...
    returnOnInvestment(investmentPrice, rewardPrice_year,...
    oneOnNPortfolio);
    %%%%%%%%%%%%%%%%%%%%%%% 1/N portfolio variance %%%%%%%%%%%%%%%%%%%%%%%%
    oneOnNVar(dataSetNumber)=...
                        calculatePortfolioVarience(oneOnNPortfolio, Sigma);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    frontierPlotTest=...
          [betaMin<ExpectedOneOnNReturns max(muvec)>ExpectedOneOnNReturns];
    if frontierPlotTest(1, :)==0
        betaMin=ExpectedOneOnNReturns;
    end
    beta=ExpectedOneOnNReturns;
    betavec=beta*ones(numPrices-1, 1);
    
    
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%% calculate unconstrained Markowitz portfolios %%%%%%%%%%
% %     MarkowitzPortfolio=ComputePortfolio(M, Sigma, beta); 
% %     
% %     %Do NOT need to calculate ANY OTHER expected returns! By design, any 
% %     %portfolio calculated from this point out is supposed to have an
% %     %expected return of beta!!!!!
% % 
% %     %%%%%%%%%%%%%%%%%% Markowitz portfolio variance %%%%%%%%%%%%%%%%%%%%%%%
% %     MarkowitzVar(dataSetNumber)=...
% %                     calculatePortfolioVarience(MarkowitzPortfolio, Sigma);
% %     %%%%%%%%%%%%%%%% calculate returns on Markowitz portfolios %%%%%%%%%%%%
% %     realMarkowitzReturns_month(dataSetNumber)=...
% %     returnOnInvestment(investmentPrice, rewardPrice_1month,...
% %     MarkowitzPortfolio);
% %     realMarkowitzReturns_3months(dataSetNumber)=...
% %     returnOnInvestment(investmentPrice, rewardPrice_3months,...
% %     MarkowitzPortfolio);
% %     realMarkowitzReturns_year(dataSetNumber)=...
% %     returnOnInvestment(investmentPrice, rewardPrice_year,...
% %     MarkowitzPortfolio);
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %               L1 REGULARIZATION GOES HERE!
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%% Backwards deletion strategy %%%%%%%%%%%%%%%%%%%%%%
%     sparsity_number=11;
%     for desiredSparsity=100:-10:10
%         sparsity_number=sparsity_number-1;
%         if sparsity_number==10
%             nonSupport=[];
%         end
%         backwardsPort=abs(MarkowitzPortfolio);          
%         while nnz(backwardsPort)>desiredSparsity
%             forDeletion=deleteThisIndex(backwardsPort, nonSupport);
%             nonSupport=[nonSupport, forDeletion]; %add the index to the non-support
%             M1=restrictM(nonSupport, M); %M1 is M with rows in the set ns 
%                                %deleted  
%             Q1=restrictSigma(Sigma, nonSupport); %Q1 is sigma with rows and 
%                                        %cols ns deleted
%             y=ComputePortfolio(M1, Q1, beta);%y is the portfolio 
%                                            %obtained from ignoring
%                                            %the contribution of 
%                                            %assets in ns
%             backwardsPort=reInsert_ns(y, nonSupport);%insert 0s into the non-support 
%                                        %elements of y  
%         end
%         selectedAssets=find(backwardsPort);
%         backwardsSelectedAssets{dataSetNumber, 11-sparsity_number}=...
%             names(selectedAssets);
%         %%%%%%%%%%%%%%%%%% backwards portfolio variance %%%%%%%%%%%%%%%%%%%%
%         BackwardsVar(dataSetNumber, 11-sparsity_number)=...
%                     calculatePortfolioVarience(backwardsPort, Sigma);     
%         %%%%%%%%%%%%%%%%%% backwards portfolio returns %%%%%%%%%%%%%%%%%%%%
%         realBackwardsReturns_month(dataSetNumber, 11-sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_1month,...
%         backwardsPort);
%         realBackwardsReturns_3months(dataSetNumber, 11-sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_3months,...
%         backwardsPort);
%         realBackwardsReturns_year(dataSetNumber, 11-sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_year,...
%         backwardsPort);          
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    sparsity_number=0;
    for desiredSparsity=10:10:100
        sparsity_number=sparsity_number+1
        
        %%%%%%%%%%%%%%%%%%%% MIQP solution by YALMIP %%%%%%%%%%%%%%%%%%%%%%
        exactSolution=MIP(Sigma, beta, muvec, desiredSparsity)
        selectedAssets=find(exactSolution);
        ExactSelectedAssets{dataSetNumber, sparsity_number}=...
            names(selectedAssets);
        %%%%%%%%%%%%%%%%%% exact portfolio variance %%%%%%%%%%%%%%%%%%%%
        ExactVar(dataSetNumber, sparsity_number)=...
                    calculatePortfolioVarience(exactSolution, Sigma);     
        %%%%%%%%%%%%%%%%%% exact portfolio returns %%%%%%%%%%%%%%%%%%%%
        ExactReturns_month(dataSetNumber, sparsity_number)=...
        returnOnInvestment(investmentPrice, rewardPrice_1month,...
        exactSolution);
        ExactReturns_3months(dataSetNumber, sparsity_number)=...
        returnOnInvestment(investmentPrice, rewardPrice_3months,...
        exactSolution);
        exactReturns_year(dataSetNumber, sparsity_number)=...
        returnOnInvestment(investmentPrice, rewardPrice_year,...
        exactSolution);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%% Forwards selection strategy %%%%%%%%%%%%%%%%%%
%         if desiredSparsity==10
%             S=getBestIndex(muvec);
%         end
%         while nnz(S)<desiredSparsity
%             nnz(S);
%             S=listInds(S, numAssets); %want to update this until 
%                                   %numrows(allInds)=desiredSparsity
%             gmvVec=getGMVreturns(Sigma, S, muvec);
%             [~, bestInd]=max(gmvVec);
%             S=S(bestInd, :);         
%             M2=M(S, :);
%             Q2=expandSigma(S, Sigma);
%             y=ComputePortfolio(M2, Q2, beta);
%             ForwardsPort=buildForwardsX(y,S, numAssets);
%         end
%         selectedAssets=find(ForwardsPort);
%         ForwardsSelectedAssets{dataSetNumber, sparsity_number}=...
%             names(selectedAssets);
%         %%%%%%%%%%%%%%%%%% forwards portfolio variance %%%%%%%%%%%%%%%%%%%%
%         ForwardsVar(dataSetNumber, sparsity_number)=...
%                     calculatePortfolioVarience(ForwardsPort, Sigma);     
%         %%%%%%%%%%%%%%%%%% forwards portfolio returns %%%%%%%%%%%%%%%%%%%%
%         realForwardsReturns_month(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_1month,...
%         ForwardsPort);
%         realForwardsReturns_3months(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_3months,...
%         ForwardsPort);
%         realForwardsReturns_year(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_year,...
%         ForwardsPort);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%% Pursuit selection setup%%%%%%%%%%%%%%%%%%%
%         A=retmat;
%         W=sqrt(diag(A'*A));
%         for k=1:1:numAssets, 
%             A(:,k)=A(:,k)/W(k); 
%         end;
%         b=beta*ones(numPrices-1, 1);
% 
%         %%%%%%%%%%%%%%%%%%%% all selection strategy %%%%%%%%%%%%%%%%%%%
%         %idea:            
%         %for each dataset
%         %    for each target return
%         %        for each desired sparsity
%         %            use PA to select a K-sparse approximant to the vector 
%         %            b(=beta) with respect to the dictionary formed by the 
%         %            matrix of all returns(=retmat)
%         %        end
%         %        use formula to reclaibrate the K-sparse approximant into
%         %        a portfolio with target return beta
%         %    end
%         %    calculate the returns over the time horizons
%         %end
%         %%%%%%%%%%%%%%%%%%%%% Implement Part (1) %%%%%%%%%%%%%%%%%%%%%%
%         OMPX=zeros(numAssets, 1);
%         r=b;    
%         if desiredSparsity==10
%             SS=[];%mark-about 3
%         end
%         while nnz(OMPX)<desiredSparsity
%             if  r'*r>threshold,
%                 Z=abs(A'*r);
%                 posZ=find(Z==max(Z));
%                 SS=sort([SS,posZ(1)]);
%                 r=b-A(:,SS)*pinv(A(:,SS))*b;    
%             end;
%             OMPX(SS)=pinv(A(:,SS))*b;
%         end
%         %%%%%%%%%%%%%%%%%%%%% Implement Part (2) %%%%%%%%%%%%%%%%%%%%%%
%         M3=M(SS, :);
%         Q3=expandSigma(SS, Sigma);
%         y=ComputePortfolio(M3, Q3, beta);
%         OMPX(SS)=y;
%         selectedAssets=find(OMPX);
%         OMPSelectedAssets{dataSetNumber, sparsity_number}=...
%             names(selectedAssets);
%         %%%%%%%%%%%%%%%%%%%%% Implement Part (3) %%%%%%%%%%%%%%%%%%%%%%
%         OMPVar(dataSetNumber, sparsity_number)=...
%                 calculatePortfolioVarience(OMPX, Sigma);     
%         %%%%%%%%%%%%%%%%%% OMP portfolio returns %%%%%%%%%%%%%%%%%%%%
%         realOMPReturns_month(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_1month,...
%         OMPX);
%         realOMPReturns_3months(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_3months,...
%         OMPX);
%         realOMPReturns_year(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_year,...
%         OMPX);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         
%         %%%%%%%%%%%%%%%%%%%%%%%% MP selection strategy %%%%%%%%%%%%%%%%%%%%%
%         MPX=zeros(numAssets,1);
%         r=b;
%         while nnz(MPX)<desiredSparsity
%             if r'*r>threshold,
%                 Z=abs(A'*r);
%                 posZ=find(Z==max(Z),1);
%                 MPX(posZ)=MPX(posZ)+A(:,posZ)'*r;
%                 r=r-A(:,posZ)*A(:,posZ)'*r;
%             end;
%             SS=find(abs(MPX)>1e-8)';
%         end
%         M4=M(SS, :);
%         Q4=expandSigma(SS, Sigma);
%         y=ComputePortfolio(M4, Q4, beta);
%         MPX(SS)=y;   
%         selectedAssets=find(MPX);
%         MPVar(dataSetNumber, sparsity_number)=...
%                 calculatePortfolioVarience(MPX, Sigma);
%         MPSelectedAssets{dataSetNumber, sparsity_number}=...
%             names(selectedAssets);     
%         %%%%%%%%%%%%%%%%%% MP portfolio returns %%%%%%%%%%%%%%%%%%%%
%         realMPReturns_month(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_1month,...
%         MPX);
%         realMPReturns_3months(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_3months,...
%         MPX);
%         realMPReturns_year(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_year,...
%         MPX);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              
%         %%%%%%%%%%%%%%%%%%%%%% LSOMP selection strategy %%%%%%%%%%%%%%%%%%%
%         
%         r=b;
%         %if sparsity_number==1 
%             SS=[];
%             %nnz(SS)
%         %end
%         while nnz(SS)<desiredSparsity
%             if r'*r>threshold,
%                 Z=zeros(numAssets,1);
%                 for jj=1:1:numAssets
%                     SStemp=[SS,jj];
%                     rtemp=b-A(:,SStemp)*pinv(A(:,SStemp))*b;
%                     Z(jj)=rtemp'*rtemp;
%                 end;
%                 posZ=find(Z==min(Z),1);
%                 SS=sort([SS,posZ(1)]);
%                 r=b-A(:,SS)*pinv(A(:,SS))*b;    
%             end;
%             LSOMPX=zeros(numAssets, 1);
%             LSOMPX(SS)=pinv(A(:,SS))*b;
%         end
%         M5=M(SS, :);
%         Q5=expandSigma(SS, Sigma);
%         y=ComputePortfolio(M5, Q5, beta);
%         LSOMPX(SS)=y; 
%         LSOMPVar(dataSetNumber, sparsity_number)=...
%             calculatePortfolioVarience(LSOMPX, Sigma);
%         selectedAssets=find(LSOMPX);
%         LSOMPSelectedAssets{dataSetNumber, sparsity_number}=...
%             names(selectedAssets)   ;
%         %%%%%%%%%%%%%%%%%%%% LSOMP portfolio returns %%%%%%%%%%%%%%%%%%%%%%
%         realLSOMPReturns_month(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_1month,...
%         LSOMPX);
%         realLSOMPReturns_3months(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_3months,...
%         LSOMPX);
%         realLSOMPReturns_year(dataSetNumber, sparsity_number)=...
%         returnOnInvestment(investmentPrice, rewardPrice_year,...
%         LSOMPX);
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        
% 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%             %%%%%%%%%%%%%%%%% tresholding selection strategy %%%%%%%%%%%%%%
%             thrX=zeros(numAssets, 1); 
%             Z=A'*b; 
%             [Za,posZ]=sort(abs(Z),'descend');
%             r=b;
%             SS=[];
%             while nnz(thrX)<desiredSparsity
%                 if r'*r>threshold, 
%                     SS=[SS,posZ(length(SS)+1)];
%                     thrX(SS)=pinv(A(:,SS))*b;
%                     r=b-A(:,SS)*thrX(SS);
%                 end;                
%             end
%             M6=M(SS, :);
%             Q6=expandSigma(SS, Sigma);
%             y=ComputePortfolio(M6, Q6, beta);
%             thrX(SS)=y;
%             thrVar(dataSetNumber, sparsity_number)=...
%             calculatePortfolioVarience(thrX, Sigma);
%             selectedAssets=find(thrX);
%             thrSelectedAssets{dataSetNumber, sparsity_number}=...
%                 names(selectedAssets)   ;
%             %%%%%%%%%%%%%%%%%%%% thr portfolio returns %%%%%%%%%%%%%%%%%%%%%%
%             realthrReturns_month(dataSetNumber, sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_1month,...
%             thrX);
%             realthrReturns_3months(dataSetNumber, sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_3months,...
%             thrX);
%             realthrReturns_year(dataSetNumber, sparsity_number)=...
%             returnOnInvestment(investmentPrice, rewardPrice_year,...
%             thrX);           
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

    

