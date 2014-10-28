function [smallDataset]= Generatedata(smallrow, smallcol)
smallDataset=randn(smallrow, smallcol);


    




% A=randn(smallrow, smallcol);
% j=randsample(1: smallrow, randi([0, smallrow]))
% if nnz(j)~=0
%     for i=1:length(j)
%         k=randsample(smallrow:smallcol, 1)
%         A(:, j(i))=rand(1)*A(:, k)
%     end
% end

    



k=randi([1 smallcol]);
A=zeros(smallrow, smallcol);
for i = 1:k, A = A + randn(smallrow,1) * randn(1,smallcol); end
smallDataset=A;