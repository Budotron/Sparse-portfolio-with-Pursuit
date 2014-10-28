 [~, initial]=max(muvec); %initialize the support with the index of the 
                                     % asset with the largest expected
                                     % return
            S=[initial]; %initial support is the index of the largest asset weight
            x=MarkowitzPortfolios{set, count};
            y=zeros(length(x), 1); 
            while nnz(y)<desiredSparsity
                for i=1:desiredSparsity
                    candidateInds=getInds(i, S, x);
                    [p, q]=size(candidateInds);
                    bestReturn=0;
                    for k=1:p %for every row of the candidates
                        tempS=candidateInds(k, :);
                        %M2=expandM(tempS, M)
                        M2=M(tempS, :);
                        Q2=expandSigma(tempS, Sigma);
                        tempy=ComputePortfolio(M2, Q2, j);
                        [a, b]=size(retmat);
                        histrets_for_testing=retmat(a-29, :);
                        tempReturn=histrets_for_testing(tempS)*tempy;
                        if tempReturn>bestReturn
                            bestReturn=tempReturn;
                            bestInds=tempS;
                        end
                    end
                    S=unique([S, bestInds])
                    y(bestInds)=tempy;
                end
            end
            expectedForwardReturns(set, count)=...
                muvec'*y
            
            
            %%%%secondtry%%%
S=[]; %initially empty support
           [~, I] = max(muvec) %initial index for support set
           S=[S, I]
           weightsTank=0; %assuming that larger weights have a greater 
                         %impact on the portfolio, we want, in each
                         %iteration to select the portfolio whose absolute
                         %value sum is the largest, and hold that index set for
                         %the next iteration. These sums are stored in the
                         %weights tank for comparison to the next portfolio
           x=zeros(length(MarkowitzPortfolios{set, count}), 1);
           while nnz(x)<desiredSparsity 
               [candidateSupports]= getIndex(S, x);%list all possible 
                        %expansions by one of the support 
               for k=1: length(candidateSupports) 
                   if k<length(candidateSupports)
                       k
                       testSupport=candidateSupports(k, :);
                       M2=M(testSupport, :); %restrict M to the support set
                       Q2=expandSigma(testSupport, Sigma);%restrict Sigma to the 
                                                 %support set
                       tempPort=ComputePortfolio(M2, Q2, j); %calculate a portfolio 
                                            %using the support set
                       [a, b]=size(retmat);
                       histrets_for_testing=retmat(a-29, :);
                       tempReturns=histrets_for_testing(testSupport)*tempPort;
                       if tempReturns>weightsTank
                           weightsTank=tempReturns
                           indexTank=testSupport;
                           portfolioTank=tempPort;
                       end
                   else
                       break
                   end
               end
           S=unique([S, indexTank] ) 
           pause
           x(S)=portfolioTank;
           end
          end