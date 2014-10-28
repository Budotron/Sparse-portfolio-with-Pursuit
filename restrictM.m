function [M1]=restrictM (ns, M)
tempM=M;
tempM(ns, :)=[];
M1=tempM;
