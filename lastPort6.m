%something to work on: how to hot-start the backwards deletion/forwards
%selection method in each iteration for each dataset

%the following code computes portfolios using 9 different methods. each
%method is applied to ten dataset, each with ten different target returns
%(though the ten returns are the same for each dataset). 

%1. set everything up. load names of assets and ten sets of daily prices
%for all assets. 

clc; clear all
warning('off', 'all');
warning;
load names %the ith name in the list corresponds to the ith column in 
           %each dataset

           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% load datasets (tempdata is rendomly generated)%%%%%

load datasets %10 sets of 300 daily proces for 440 S&P listed companies

% tempData=cell(10, 1);
% for i=1:length(tempData)
%     tempData{i}=Generatedata(70, 110);
% end
% datasets=tempData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[numPeriods, ~]=size(datasets); %period will be used to select each 
                            %dataset

                            
%2. loop through the sets of data, and compute portfolios for different 
%levels of target return for each set. All datasets have the same number of
%assets (numAssets), so the 1/N portfolio is the same in every case. each 
%dataset (period) will give one realized return value for the 1/N
%portfolio, each of which is stored in a 1 x 10 matrix (*). There are three
%return periods per dataset. We store each 1 x 10 realized return matrix 
%in a 10 x 3
%cell, the ith column of which represents the ith return period (**). 
%There is one unrestricted Markowitz portfolio per target return per
%dataset (one hundred total). These are stored in a 10 x 10 cell(?), where 
%column i of the cell holds portfolios for dataset i, and row j represents
%target return j. For each portfolio, there is one realized return per
%target return per dataset. One hundred portfolios thus generate 300 target
%returns. These results are stored in three 10 x 10 matrices(?), where matrix
%k represents the kth time horizon. The first portfolio (for dataset1,
%target retun 1) will have a realized return over the three time horizons
%stored in the upper left corner of each. Realized returns for 
%portfolio (2, 1) (for dataset1, %target return2) will be stored in the 
%(2, 1) position of each matrix (??). 
%Each restricted cardinality problem has a target cardinality per target
%return per dataset. There are 10 cardinalities per target return, and ten
%target returns per dataset. In each dataset, we claculate 10 portfolios
%(one for each cardinality) for 10 target returns (100 per dataset; 1000
%altogether). Approach: loop through datasets; in dataset i, loop through
%target returns; in target return j, loop through cardinalities. For the
%first target return, generate 10 portfolios corresponding to the 10
%cardinalities. Store each portfolio as a column in a numAssets x 10
%matrix, where column i is stores a portfolio of cardinality i (+). There
%are ten such matrices (one for each target return). Store these in a 1 x
%10 cell. Cell i corresponds to target return i. Finally, store each 1 x 10
%cell in a 10 x 1 cell, each of which is a different dataset. All porfolios
%are thus stored in a 10 x 10 cell, each row of which is a different
%dataset, and each column of which is a different target return (++). 
%Each portfolio has 3 time horizons(?), with cell i denoting horizon i.
%Each cell contains a 10 x 10 matrix (??), the (1,1) position of each of 
%which is the %realised return in period i of the portfolio of sparsity 1 
%and target return 1, and the (1, 2) position of each of which is the return
%of portfolio 1 and target return 2. 


realOneOnNReturns_day=zeros(1, numPeriods); %(*)
realOneOnNReturns_week=zeros(1, numPeriods);%(*)
realOneOnNReturns_month=zeros(1, numPeriods);%(*)
realOneOnNReturns_all=cell(1, 3);%(**)
oneOnNCalctime=zeros(1,1);
MarkowitzPortfoliosMat=cell(numPeriods, numPeriods); %(?)
realMarkowitzReturns_Day=zeros(numPeriods, numPeriods);%(?)
realMarkowitzReturns_Week=zeros(numPeriods, numPeriods);%(?)
realMarkowitzReturns_Month=zeros(numPeriods, numPeriods);%(?)
MarkowitzReturns_all=cell(1, 3); %(??) 
MarkowitzCalctimes=zeros(numPeriods, numPeriods);
MarkowitzVarMat=zeros(numPeriods, numPeriods);
BackwardsNamesMat=cell(numPeriods, numPeriods);
realBackwardsReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealBackwardReturns_Day=cell(1, numPeriods);
BackwardsNamesMat_all=cell(1, numPeriods);
realBackwardsReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealBackwardReturns_Week=cell(1, numPeriods);
realBackwardsReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealBackwardReturns_Month=cell(1, numPeriods);
BackwardsCalctimesMat=zeros(numPeriods, numPeriods);
BackwardsCalctimes=cell(1, numPeriods);
allBackwardsCalctimes=cell(1, numPeriods);
BackwardsVarMat=zeros(numPeriods, numPeriods);
BackwardsVar_all=cell(1, numPeriods);
realForwardsReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealForwardsReturns_Day=cell(1, numPeriods);
realForwardsReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealForwardsReturns_Week=cell(1, numPeriods);
realForwardsReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealForwardsReturns_Month=cell(1, numPeriods);
ForwardsNamesMat=cell(numPeriods, numPeriods);
ForwardsNamesMat_all=cell(1, numPeriods);
ForwardsVarMat=zeros(numPeriods, numPeriods);
ForwardsVar_all=cell(1, numPeriods);
realOMPReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealOMPReturns_Day=cell(1, numPeriods);
realOMPReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealOMPReturns_Week=cell(1, numPeriods);
realOMPReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealOMPReturns_Month=cell(1, numPeriods);
OMPNamesMat=cell(numPeriods, numPeriods);
OMPNamesMat_all=cell(1, numPeriods);
OMPVarMat=zeros(numPeriods, numPeriods);
OMPVar_all=cell(1, numPeriods);
realMPReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealMPReturns_Day=cell(1, numPeriods);
realMPReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealMPReturns_Week=cell(1, numPeriods);
realMPReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealMPReturns_Month=cell(1, numPeriods);
MPNamesMat=cell(numPeriods, numPeriods);
MPNamesMat_all=cell(1, numPeriods);
MPVarMat=zeros(numPeriods, numPeriods);
MPVar_all=cell(1, numPeriods);
realLSOMPReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealLSOMPReturns_Day=cell(1, numPeriods);
realLSOMPReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealLSOMPReturns_Week=cell(1, numPeriods);
realLSOMPReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealLSOMPReturns_Month=cell(1, numPeriods);
LSOMPNamesMat=cell(numPeriods, numPeriods);
LSOMPNamesMat_all=cell(1, numPeriods);
LSOMPVarMat=zeros(numPeriods, numPeriods);
LSOMPVar_all=cell(1, numPeriods);
realThrReturns_Day=zeros(numPeriods, numPeriods);%(?)
allRealThrReturns_Day=cell(1, numPeriods);
realThrReturns_Week=zeros(numPeriods, numPeriods);%(?)
allRealThrReturns_Week=cell(1, numPeriods);
realThrReturns_Month=zeros(numPeriods, numPeriods);%(?)
allRealThrReturns_Month=cell(1, numPeriods);
thrNamesMat=cell(numPeriods, numPeriods);
thrNamesMat_all=cell(1, numPeriods);
thrVarMat=zeros(numPeriods, numPeriods);
thrVar_all=cell(1, numPeriods);
threshold=1e-4;
%3. begin loop

for dataSetNumber=1:numPeriods   
    dailyPrices=datasets{dataSetNumber}; %each dataset is prepared for analysis
    [numDays,prices]=size(dailyPrices);
    lastTenth=1/10*numDays;
    purchaseDayPrices=dailyPrices(numDays-lastTenth, :);
    nextDayPrices=dailyPrices(numDays-(lastTenth-1), :);
    nextWeekPrices=dailyPrices(numDays-(lastTenth-7), :);
    nextMonthPrices=dailyPrices(numDays, :);
    retmat=makeReturnsMatrix(dailyPrices(1:numDays-lastTenth, :)); 
                                                            %R_{t+1}/R{t}-1
    [retrows retcols]=size(retmat);
    [ReturnDays, dailyReturn]=size(retmat);
    muvec=mean(retmat)'; %vector of expected returns
    Sigma=cov(retmat);   %covariance matrix of return vectors
    onesvec=ones(length(muvec), 1);
    betaMax=0.1; %largest target return
    gvx=globalMinimumVariancePortfolio(Sigma, onesvec); %portfolio from 
                                                        %ignoring
                                                        %target constraint
    betaMin=muvec'*gvx; %expected return of GMV portfolio. Any portfolio 
                        %with a higer expected return will be non-dominated
    betarange=sort([betaMin, betaMax], 'ascend');
    M=[muvec onesvec]; %concatenation for use in formula 
    backwardPortMat=zeros(prices, prices);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%   set up 1/N portfolio    %%%%%%%%%%%%%%%%%%%%%%%
    oneOnNPortfolio=repmat(1/dailyReturn, [dailyReturn,1]);   
    %%%%%%%%%%%%%%%% calculate returns on 1/N ports per period %%%%%%%%%%%%
    realOneOnNReturns_day(dataSetNumber)=...
    returnOnInvestment(purchaseDayPrices, nextDayPrices,...
    oneOnNPortfolio); 
    realOneOnNReturns_week(dataSetNumber)=...
    returnOnInvestment(purchaseDayPrices, nextWeekPrices,...
    oneOnNPortfolio);
    realOneOnNReturns_month(dataSetNumber)=...
    returnOnInvestment(purchaseDayPrices, nextMonthPrices,...
    oneOnNPortfolio);
    %%%%%%%%%%%%%%%%%%%%%%%%%% portfolio variance %%%%%%%%%%%%%%%%%%%%%%%%%
    oneOnNvar(dataSetNumber)=...
                        calculatePortfolioVarience(oneOnNPortfolio, Sigma)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    expectedReturnsVec=[min(betarange):...
            (max(betarange)-min(betarange))/9:max(betarange)];
        
    targetReturnNumber=0;
    
    %%%%%%%%%%%%%%%%%%%%% for each target return %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for beta=min(betarange):(max(betarange)-min(betarange))/9:max(betarange)
        targetReturnNumber=targetReturnNumber+1;   
        

        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%% calculate unconstrained Markowitz portfolios %%%%%%%%
        MarkowitzPortfoliosMat{dataSetNumber, targetReturnNumber}=...
            ComputePortfolio(M, Sigma, beta); %rows are datasets
                                              %columns are returns
        MarkowitzVarMat(dataSetNumber, targetReturnNumber)=...
        calculatePortfolioVarience(MarkowitzPortfoliosMat{dataSetNumber,...
        targetReturnNumber}, Sigma);
        %%%%%%%%%%%%% calculate unconstrained Markowitz returns %%%%%%%%%%%
        realMarkowitzReturns_Day(dataSetNumber, targetReturnNumber)=...
            returnOnInvestment(purchaseDayPrices, nextDayPrices,...
            MarkowitzPortfoliosMat{dataSetNumber, targetReturnNumber});
        realMarkowitzReturns_Week(dataSetNumber, targetReturnNumber)=...
            returnOnInvestment(purchaseDayPrices, nextWeekPrices,...
            MarkowitzPortfoliosMat{dataSetNumber, targetReturnNumber});
        realMarkowitzReturns_Month(dataSetNumber, targetReturnNumber)=...
            returnOnInvestment(purchaseDayPrices, nextMonthPrices,...
            MarkowitzPortfoliosMat{dataSetNumber, targetReturnNumber});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%% loop through cardinality constraints %%%%%
        sparsity_set=0;
        for desiredSparsity=10:10:100
            sparsity_set=sparsity_set+1;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%% Backwards deletion strategy %%%%%%%%%%%%%%%%

            x=abs(MarkowitzPortfoliosMat{dataSetNumber, targetReturnNumber});
            ns=[];%the non-support set, initially empty, is the set of 
                  %indices that will be set to 0
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
                  x=reInsert_ns(y, ns);%insert 0s into the non-support 
                                       %elements of y  
                  
            end
%            BackwardsCalctimesMat(targetReturnNumber, sparsity_set)=toc;
            BackwardsCalctimes{targetReturnNumber, sparsity_set}=...
                                                        BackwardsCalctimesMat;
            realBackwardsReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, x);
            realBackwardsReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, x);
            realBackwardsReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, x);
                %rows are expected values, columns are sparsity targets
            backwardsIndices=find(x);
            BackwardsNamesMat{targetReturnNumber, sparsity_set}=...
                names(backwardsIndices);
            BackwardsVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (x, Sigma)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%% Forwards selection strategy %%%%%%%%%%%%%%%%
            % idea: 
            % for each dataset
            %     for each target return
            %         for each desired sparsity
            %            (1) find the best perfoming global minimum variance 
            %            portfolio \leq
            %             desired sparsity
            %         end
            %        (2) use formula to recalibrate weights so that the 
            %        return on the portfolio matches target return
            %     end
            %     (3)calculate the return over the time horizons
            % end            
            %%%%%%%%%%%%%%%%%%%%% Implement Part (1) %%%%%%%%%%%%%%%%%%%%%%
            S=getBestIndex(muvec);
            indexTank=[];
            while nnz(S)<desiredSparsity
                nnz(S)
                S=listInds(S, prices); %want to update this until 
                                      %numrows(allInds)=desiredSparsity
                gmvVec=getGMVreturns(Sigma, S, muvec);
                [~, bestInd]=max(gmvVec);
                S=S(bestInd, :);
            end
            solutionSparsity=nnz(unique(S));
            [desiredSparsity solutionSparsity]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%% Implement Part (2) %%%%%%%%%%%%%%%%%%%%%%
            M2=M(S, :);
            Q2=expandSigma(S, Sigma);
            y=ComputePortfolio(M2, Q2, beta);
            x=buildForwardsX(y,S, prices);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%% Implement Part (3) %%%%%%%%%%%%%%%%%%%%%%
            realForwardsReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, x)
            realForwardsReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, x);
            realForwardsReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, x);
            ForwardsIndices=find(x);
            ForwardsNamesMat{targetReturnNumber, sparsity_set}=...
                 names(ForwardsIndices)
            ForwardsVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (x, Sigma)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Pursuit selection setup%%%%%%%%%%%%%%%%%%%
            A=retmat;
            W=sqrt(diag(A'*A));
            for k=1:1:retcols, 
                A(:,k)=A(:,k)/W(k); 
            end;
            b=beta*ones(retrows, 1);
            
            %%%%%%%%%%%%%%%%%%%% OMP selection strategy %%%%%%%%%%%%%%%%%%%
            %idea:            
            %for each dataset
            %    for each target return
            %        for each desired sparsity
            %            use OMP to select a K-sparse approximant to the vector 
            %            b(=beta) with respect to the dictionary formed by the 
            %            matrix of all returns(=retmat)
            %        end
            %        use formula to reclaibrate the K-sparse approximant into
            %        a portfolio with target return beta
            %    end
            %    calculate the returns over the time horizons
            %end
            %%%%%%%%%%%%%%%%%%%%% Implement Part (1) %%%%%%%%%%%%%%%%%%%%%%
            OMPX=zeros(prices, 1);
            r=b;
            SS=[];
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
            %%%%%%%%%%%%%%%%%%%%% Implement Part (3) %%%%%%%%%%%%%%%%%%%%%%
            realOMPReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, OMPX);
            realOMPReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, OMPX);
            realOMPReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, OMPX);
            OMPIndices=find(OMPX);
            OMPNamesMat{targetReturnNumber, sparsity_set}=...
                 names(OMPIndices);
            OMPVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (OMPX, Sigma);           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%% MP selection strategy %%%%%%%%%%%%%%%%%%%
            MPX=zeros(prices,1);
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
            realMPReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, MPX);
            realMPReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, MPX);
            realMPReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, MPX);
            MPIndices=find(MPX);
            MPNamesMat{targetReturnNumber, sparsity_set}=...
                 names(MPIndices);
            MPVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (MPX, Sigma);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             
            %%%%%%%%%%%%%%%%%%% LSOMP selection strategy %%%%%%%%%%%%%%%%%%
            LSOMPX=zeros(prices, 1);
            r=b;
            SS=[];
            while nnz(LSOMPX)<desiredSparsity
                if r'*r>threshold,
                    Z=zeros(retcols,1);
                    for jj=1:1:retcols
                        SStemp=[SS,jj]; 
                        rtemp=b-A(:,SStemp)*pinv(A(:,SStemp))*b;
                        Z(jj)=rtemp'*rtemp;
                    end;
                    posZ=find(Z==min(Z),1);
                    SS=sort([SS,posZ(1)]);
                    r=b-A(:,SS)*pinv(A(:,SS))*b;    
                end;
                LSOMPX(SS)=pinv(A(:,SS))*b;
            end
            nnz(LSOMPX);
            M5=M(SS, :);
            Q5=expandSigma(SS, Sigma);
            y=ComputePortfolio(M5, Q5, beta);
            LSOMPX(SS)=y;
            realLSOMPReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, LSOMPX);
            realLSOMPReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, LSOMPX);
            realLSOMPReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, LSOMPX);
            LSOMPIndices=find(LSOMPX);
            LSOMPNamesMat{targetReturnNumber, sparsity_set}=...
                 names(LSOMPIndices);
            LSOMPVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (LSOMPX, Sigma)           

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%% tresholding selection strategy %%%%%%%%%%%%%%
            thrX=zeros(prices, 1); 
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
            realThrReturns_Day(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextDayPrices, thrX);
            realThrReturns_Week(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextWeekPrices, thrX);
            realThrReturns_Month(targetReturnNumber, sparsity_set)=...
                returnOnInvestment(purchaseDayPrices, nextMonthPrices, thrX);
            thrIndices=find(thrX);
            thrNamesMat{targetReturnNumber, sparsity_set}=...
                 names(thrIndices);
            thrVarMat(targetReturnNumber,sparsity_set)=...
                                    calculatePortfolioVarience (thrX, Sigma)           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         end
     end
    %%%%%%%%%%%%%%%%%%%%%% store all 1/N returns per period%%%%%%%%%%%%%%%
    realOneOnNReturns_all{1}=realOneOnNReturns_day;
    realOneOnNReturns_all{2}=realOneOnNReturns_week;
    realOneOnNReturns_all{3}=realOneOnNReturns_month;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%% store all Markowitz returns per period %%%%%%%%%%%%%
    realMarkowitzReturns_all{1}=realMarkowitzReturns_Day
    realMarkowitzReturns_all{2}=realMarkowitzReturns_Week
    realMarkowitzReturns_all{3}=realMarkowitzReturns_Month
    allBackwardsCalctimes{dataSetNumber}=BackwardsCalctimes
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% store all backwards expected returns per %%%%%%%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealBackwardReturns_Day{dataSetNumber}=realBackwardsReturns_Day
    allRealBackwardReturns_Week{dataSetNumber}=realBackwardsReturns_Week
    allRealBackwardReturns_Month{dataSetNumber}=realBackwardsReturns_Month
    BackwardsNamesMat_all{dataSetNumber}=BackwardsNamesMat
    allBackwardsCalctimes{dataSetNumber}=BackwardsCalctimes
    BackwardsVar_all{dataSetNumber}=BackwardsVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%% store all forwards expected returns per %%%%%%%%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealForwardsReturns_Day{dataSetNumber}=realForwardsReturns_Day
    allRealForwardsReturns_Week{dataSetNumber}=realForwardsReturns_Week
    allRealForwardsReturns_Month{dataSetNumber}=realForwardsReturns_Month
    ForwardsNamesMat_all{dataSetNumber}=ForwardsNamesMat
    ForwardsVar_all{dataSetNumber}=ForwardsVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% store all OMP expected returns per %%%%%%%%%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealOMPReturns_Day{dataSetNumber}=realOMPReturns_Day
    allRealOMPReturns_Week{dataSetNumber}=realOMPReturns_Week
    allRealOMPReturns_Month{dataSetNumber}=realOMPReturns_Month
    OMPNamesMat_all{dataSetNumber}=OMPNamesMat
    OMPVar_all{dataSetNumber}=OMPVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%% store all MP expected returns per %%%%%%%%%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealMPReturns_Day{dataSetNumber}=realMPReturns_Day
    allRealMPReturns_Week{dataSetNumber}=realMPReturns_Week
    allRealMPReturns_Month{dataSetNumber}=realMPReturns_Month
    MPNamesMat_all{dataSetNumber}=MPNamesMat
    MPVar_all{dataSetNumber}=MPVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%% store all LSOMP expected returns per %%%%%%%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealLSOMPReturns_Day{dataSetNumber}=realLSOMPReturns_Day
    allRealLSOMPReturns_Week{dataSetNumber}=realLSOMPReturns_Week
    allRealLSOMPReturns_Month{dataSetNumber}=realLSOMPReturns_Month
    LSOMPNamesMat_all{dataSetNumber}=LSOMPNamesMat
    LSOMPVar_all{dataSetNumber}=LSOMPVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%% store all thresholding expected returns per %%%%%%%%%%%%
    %%%%% cardinality in a cell row corresponding to each time horizon %%%% 
    allRealThrReturns_Day{dataSetNumber}=realThrReturns_Day
    allRealThrReturns_Week{dataSetNumber}=realThrReturns_Week
    allRealThrReturns_Month{dataSetNumber}=realThrReturns_Month
    thrNamesMat_all{dataSetNumber}=thrNamesMat
    thrVar_all{dataSetNumber}=thrVarMat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


end
save realOneOnNReturns_day
save realOneOnNReturns_week 
save realOneOnNReturns_month
save realOneOnNReturns_all
save oneOnNvar 
save realMarkowitzReturns_all
save MarkowitzVarMat
save allRealBackwardReturns_Day
save allRealBackwardReturns_Week
save allRealBackwardReturns_Month
save BackwardsNamesMat_all
save allBackwardsCalctimes
save BackwardsVar_all
save allRealForwardsReturns_Day{dataSetNumber}
save allRealForwardsReturns_Week
save allRealForwardsReturns_Month
save ForwardsNamesMat_all{dataSetNumber}
save ForwardsVar_all
save allRealOMPReturns_Day
save allRealOMPReturns_Week
save allRealOMPReturns_Month
save OMPNamesMat_all
save OMPVar_all
save allRealMPReturns_Day
save allRealMPReturns_Week
save allRealMPReturns_Month
save MPNamesMat_all
save MPVar_all
save allRealLSOMPReturns_Day
save allRealLSOMPReturns_Week
save allRealLSOMPReturns_Month
save LSOMPNamesMat_all 
save LSOMPVar_all
save allRealThrReturns_Day
save allRealThrReturns_Week
save allRealThrReturns_Month
save thrNamesMat_all
save thrVar_all



