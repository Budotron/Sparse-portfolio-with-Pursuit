sig= 5 + (80-5).*rand(100,1);
mu= -20 + (20--20).*rand(100,1);
S0= 5 +(200-5).*rand(100,1)
nsims=1
steps=1
dt=100
for i=1:length(mu)
    nu = mu(i) - sig(i)*sig(i)/2;
    S = S0(i)*[ones(1,nsims); ...
            cumprod(exp(nu*dt+sig*sqrt(dt)*randn(steps,nsims)),1)]
        A(i, :)=max(S)
        keyboard
end
