function [Q1] =restrictSigma(Sigma, ns)
tempQ=Sigma;
tempQ(ns, :)=[];
tempQ(:, ns)=[];
Q1=tempQ;
