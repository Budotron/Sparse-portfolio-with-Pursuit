function [retmat] = makeReturnsMatrix (dailyPrices)
[m,n]=size(dailyPrices);
retmat=zeros(m-1, n);
for i=1:m-1
    retmat(i, :)=(dailyPrices(i+1, :)-dailyPrices (i, :))./dailyPrices(i, :);
end