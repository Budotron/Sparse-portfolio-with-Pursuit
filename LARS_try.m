%n=30; m=50; 
Smax=15; Exper=200; sigma=1e-2;
load('datasets.mat');
A=datasets{1,1};
A=makeReturnsMatrix(A)
[n, m]=size(A)
%A=randn(n,m);
W=sqrt(diag(A'*A));
for k=1:1:m, 
    A(:,k)=A(:,k)/W(k); 
end; 
b=0.035*ones(n, 1);
Er2=zeros(Smax,Exper,3); 
ErS=zeros(Smax,Exper,2); 
Card=zeros(Smax,Exper,2);
for S=1:1:Smax,
        [xLARS,lambda]=Chapter_05_LARS(A,b,0);
        xLARS
        lambda
        xLARS=xLARS'; 
        choice=find(sqrt(mean((A*xLARS-b*ones(1,n)).^2,1))<=sqrt(2)*sigma,1);
        Sopt=find(xLARS(:,choice))'; 
        Sopt
        %ErS(S,experiment,2)=(max(S,length(Sopt))-...
                      % length(intersect(Sopt,pos(1:S))))/max(S,length(Sopt)); 
        %Er2(S,experiment,2)=mean((xLARS(:,choice)-x0).^2)/(x0'*x0); 
        Card(S,experiment,2)=length(Sopt);
end
xLARS
lambda
