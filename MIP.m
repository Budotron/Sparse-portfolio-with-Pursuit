function [exactSolution]=MIP(Sigma, beta, muvec, desiredSparsity)
[n, ~]=size(Sigma);
w = sdpvar(n,1);
mutarget = beta;
F = [sum(w) == 1,  muvec'*w == mutarget]
F = [F, nnz(w) <= desiredSparsity];
solvesdp(F,w'*Sigma*w)
exactSolution=double(w);