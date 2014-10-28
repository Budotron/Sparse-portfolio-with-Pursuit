function [BNBX]=branchAndBoundSolve(S, mutarget, mu,  card)
[n, ~]=size(S)
w = sdpvar(n,1);
F = [sum(w) == 1,  mu*w == mutarget];
card
F = [F, nnz(w) <= card];
options = sdpsettings('solver', 'bnb');
solvesdp(F,w'*S*w, options);
BNBX=double(w);
BNBX(BNBX<=1e-4)=0;
