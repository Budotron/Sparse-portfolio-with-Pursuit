clc; clear all;
load datapoints; %load 501 historical prices for 491 S&P listed companies
load names; %the ith name in the list corresponds to the ith column in datapoints
[p, q]= size(datapoints);
sample=datapoints(p-130, :); %take the 100 prices of every asset leaving 
                                  %days out of sample to test portfolio
                                  %perfomance
retmat1=makeReturnsMatrix(sample); %R_{t+1}/R{t}-1
[j, k]=size(retmat1)
retmat=retmat(1:j-30, :)
% load mu; load sigma
% mu
% sigma
% muvec=mu
% Sigma=sigma
muvec=mean(retmat)'; %vector of expected returns
Sigma=cov(retmat);   %covariance matrix of return vectors
onesvec=ones(length(muvec), 1);
[betaMax, I]=max(muvec);%the largest asset expected return is taken as the 
                        %the largest target expected portfolio return
%%%
betaMax=.1 %bc previous betaMax is too small
%%%
gvx=globalMinimumVariancePortfolio(Sigma, onesvec); %portfolio from ignoring
                                                    %target constraint
betaMin=muvec'*gvx; %expected return of GMV portfolio. Any portfolio with 
                    %a higer expected return will be non-dominated
betarange=sort([betaMin, betaMax], 'ascend');
M=[muvec onesvec]; %concatenation for use in formula
count=0;
portfolios=zeros(k, 10);
for i=min(betarange):(max(betarange)-min(betarange))/10:max(betarange)
    %the following calculates the expected returns for all portfolios
    %calculated with target expected return in the range [betaMin, betaMax]
    count=count+1;
    portfolios(:, count)=ComputePortfolio(M, Sigma, i)
    usualExpRets(count)=muvec'*portfolios(:, count); %compute the expected return for each portfolio
                             %produced
    %realdayrets(count)=datapoints(p-29, :)*portfolios(:, count);                        
    UsualPortVar=calculatePortfolioVarience (portfolios(:, count), Sigma);
    UsualSD(count)=sqrt(UsualPortVar); %compute the standdard deviation for each
                                 %portfolio produced
    
    %The following calculates the expected and real returns on portfolios 
    %the 1/N strategy
    
    oneOnN=repmat(1/k, [k,1]);
    oneOnNExpRets(count)=muvec'*oneOnN;
    realOneOnNrets(count)=datapoints(p-29, :)*oneOnN;
    
end
% hold on
% plot(stddev, exprets, '*')
% plot(stddev, realdayrets, 'r*')