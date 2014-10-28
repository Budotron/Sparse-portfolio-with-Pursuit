function [hold]=bopRepeats(S)
[m,n]=size(S);
hold=[];
for I=1:m
    check=nnz(unique(S(I, :)));
    if check<n
        hold=[hold I];
    end
end