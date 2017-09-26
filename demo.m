clear all, close all
n=5;
m=5;
p=3;
d=5;
L=1024; % L is total number of samples
k=floor((L)/(p+nchoosek(d+m-1,m-1)));
r=k*nchoosek(d+m-1,m-1)-(k-1);
N=r+k*p;    % N is number of columns of block Hankel matrices
    
% generate stable A matrix
A=rand(n);
A=A-max(abs(eig(A)))*eye(n);
A=((eye(n)+A)/(eye(n)-A))*.5;
B=randn(n,m^d);
B(:,1)=zeros(n,1);
C=rand(p,n);
D=randn(p,m^d);
D(:,1)=zeros(p,1);

u=randn(L*2,m-1);
y=simssd(A,B,C,D,zeros(n,1),u,d);
    
[Ahat,BDhat,Chat]=TNmoesp(y,u,d);

% verify system matrices
u2=randn(L*5,m-1);
y2=simssd(A,B,C,D,zeros(n,1),u2,d);
yhat=simTNss(Ahat,BDhat,Chat,zeros(size(Ahat,1),1),u2);    
for i=1:p
	figure
	plot(y2(:,i),'-og');grid on,hold on,plot(yhat(:,i))
end
norm(y2(:)-yhat(:))/norm(y2(:))


