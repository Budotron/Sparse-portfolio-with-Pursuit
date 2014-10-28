function [allInds]=listInds(S, prices)
list=[1:prices]';
allInds=[repmat(S, [length(list), 1]), list];