%return on investment
function [roi]=returnOnInvestment(purchaseDayPrices, futurePrice, x)
investment=purchaseDayPrices*x;
payoff=futurePrice*x;
roi=(payoff-investment)/investment;