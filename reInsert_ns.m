function [newx]=reInsert_ns(y, ns)
newx=zeros(length(y)+nnz(ns), 1);
supp=[1:length(newx)];
supp(ns)=[];
newx(supp)=y;
    

            
    
   