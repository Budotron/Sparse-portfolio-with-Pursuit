clc
clear all
numAssets=100
count=0
for i=1:10
    count=count+1
    load datapoints
    allData=datapoints(length(datapoints)-399:length(datapoints), :);
    window=1:300;
    numPrices=300;
    threshold=1e-8;
    thisPeriodPrices=datapoints(window, :);
    choose=randsample(440, 100);
    thisPeriodPrices=thisPeriodPrices(:, choose);
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
    oneOnNPortfolio=repmat(1/numAssets, [numAssets,1]);  
    ExpectedOneOnNReturns=muvec'*oneOnNPortfolio;
    beta=ExpectedOneOnNReturns;
    sparsecount=0;
    for desiredSparsity=10:10:30
        sparsecount=sparsecount+1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%% BNB solution by YALMIP %%%%%%%%%%%%%%%%%%%%%%
        BNBX=...
        branchAndBoundSolve(Sigma, beta, muvec', desiredSparsity);
        capture(count, sparsecount)=nnz(BNBX)
        
        assetsChosenBNBX=names(find(BNBX));
    %     sparsity_number=0; 
    %     for desiredSparsity=10:10:50 
    %         sparsity_number=sparsity_number+1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% Pursuit selection setup%%%%%%%%%%%%%%%%%%%
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
            OMPSelectedAssets=...
                names(selectedAssets);
            OMPOverlap(count, sparsecount)=sum(nnz(ismember(assetsChosenBNBX, OMPSelectedAssets)))
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
            selectedAssets=find(MPX);
            MPSelectedAssets=names(selectedAssets);     
            MPOverlap(count, sparsecount)=sum(nnz(ismember(assetsChosenBNBX, MPSelectedAssets)))
                        %%%%%%%%%%%%%%%%%%%% LSOMP selection strategy %%%%%%%%%%%%%%%%%%%

            r=b;
            %if sparsity_number==1 
                SS=[];
                nnz(SS)
            %end
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
            M5=M(SS, :);
            Q5=expandSigma(SS, Sigma);
            y=ComputePortfolio(M5, Q5, beta);
            LSOMPX(SS)=y; 
            selectedAssets=find(LSOMPX)
            LSOMPSelectedAssets=...
                names(selectedAssets)   ;
            LSOMPOverlap(count, sparsecount)=sum(nnz(ismember(assetsChosenBNBX, LSOMPSelectedAssets)))

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
            selectedAssets=find(thrX)
            thrSelectedAssets=...
                names(selectedAssets)   ;
            thrOverlap(count, sparsecount)=sum(nnz(ismember(assetsChosenBNBX, thrSelectedAssets)))

    end
 end