function [forDeletion]=deleteThisIndex(x, ns)
x(ns)=0;
[y, I]=min(abs(x(x~=0)));
P=find(abs(x)==y);
forDeletion=P;

% x=abs(x);

% [~, forDeletion]=min(x(x~=0));


