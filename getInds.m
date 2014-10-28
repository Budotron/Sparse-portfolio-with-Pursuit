function [candidateInds]=getInds (S, MarkowitzPortfoliosMat, dataSet,...
                                                        targetReturnNumber)
candidateIndsLength=length(MarkowitzPortfoliosMat{dataSet,targetReturnNumber});
allInds=[1:candidateIndsLength]';
candidateInds=[(repmat(S, [candidateIndsLength, 1])), allInds];


%below is best so far
% allInds=length(x)
% A=combinator(allInds, i, 'c');
% B=ismember(A, S);
% C=sum(B, 2);
% inds=find(C==i-1);
% candidateInds=A(IndsSet, :);