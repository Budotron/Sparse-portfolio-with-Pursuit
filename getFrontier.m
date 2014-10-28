function [efficientFrontier]=...
    getFrontier(gvx, mrx, mu, S, betaMin, betaMax)
grid=1:-.1:-1
minVar=gvx'*S*gvx;;
betamaxVar=mrx'*S*mrx;
covvar=gvx'*S*mrx
count=0;
for i=grid
    count=count+1
    effMu(count)=i*betaMin+(1-i)*betaMax
    effVar(count)=i^2*minVar+(1-i)^2*betamaxVar+2*i*(1-i)*covvar
end
vareach=(diag(S));
hold on
efficientFrontier=plot(effVar, effMu )
plot(vareach, mu, '*r')
xlabel('Variance')
ylabel('Expected Return')

    
