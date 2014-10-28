function [candidateSupports]= getIndex(S, x_nought)
listIndices=[1:length(x_nought)]';
currentSupport=repmat(S, [length(x_nought), 1]);
candidateSupports=[currentSupport listIndices];