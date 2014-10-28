function[portvar]=calculatePortfolioVarience (assetWeights, Sigma)
portvar=assetWeights'*Sigma*assetWeights;