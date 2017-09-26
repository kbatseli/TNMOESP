function [A,B,C,D]=moespd(y,u,d,varargin)
% [A,B,C,D]=moespd(y,u,d) or [A,B,C,D,e]=moespd(y,u,d,n)
% --------------------------------------------------------
% Matrix implementation of the MOESP algorithm for the identification of
% a specific polynomial state space model. 
%
% A 		= 	matrix, n x n state update matrix,
%
% B 		=	matrix, n x m^d input-to-state matrix, m is the number of inputs,
%
% C 		=	matrix, p x n state output model matrix,
%
% D 		= 	matrix, p x m^d input-to-output matrix,
%
% y         =   matrix, an L x p matrix where y(:,k) is the kth output,
%
% u 		=	matrix, L x m matrix where u(:,k) is the kth input,
%
% d 		=	scalar, maximal total degree of polynomial nonlinearity,
%
% n 		=	scalar, optional: maximal system order (=size of A).
%
% Reference
% ---------
%
% Tensor network subspace identification of specific polynomial state space models
%
% 2017, Kim Batselier, Ching Yun Ko, Ngai Wong

% correct for the affine entry in the extended u vector
u=[ones(size(u,1),1),u];
[L,m]=size(u);
p=size(y,2);

% k is determined from how many samples were measured.
% we set dim(U,2)=N=rank(U)
% k=min(floor((L+1)/(p+nchoosek(d+m-1,m-1))),150);
k=floor((L)/(p+nchoosek(d+m-1,m-1)));
r=k*nchoosek(d+m-1,m-1)-(k-1);
N=r+k*p;

% construct block Hankel matrices
Y=zeros(k*p,N);
U=zeros(k*m^d,N);
for i=1:N
    temp=u(i:i+k-1,:)';
    for j=2:d
        temp=khatr(temp,u(i:i+k-1,:)');
    end        
    U(:,i)=reshape(temp,[k*m^d,1]);
    Y(:,i)=reshape(y(i:i+k-1,:)',[k*p,1]);
end
% SVD-based, N=r+kp
if k*m^d > N
    [Uu,Su,Vu]=svd(U,'econ');
else
    [Uu,Su,Vu]=svd(U);
end
% inspect persistency of excitation of order k on input by checking the
% rank-gap
s=diag(Su);
disp(['Rank gap of U matrix is: ' num2str(s(r)/s(r+1))])
L11=Uu(:,1:r)*Su(1:r,1:r);
L21=Y*Vu(:,1:r);
L22=Y*Vu(:,r+1:end);

% % QR-based
% X=qr([U;Y]',0);
% L=triu(X)';
% L11=L(1:k*m^d,1:k*m^d);
% L21=L(k*m^d+1:k*m^d+k*p,1:k*m^d);
% L22=L(k*m^d+1:k*m^d+k*p,k*m^d+1:end);

[U,S,~]=svd(L22,'econ');
% determine system order n as numerical rank of L22
if ~isempty(varargin)
    n=varargin{1};
else
    tol=numel(L22)*eps;
    s=diag(S)/S(1,1);
    n=min(sum(s>tol),100);
end
O=U(:,1:n)*diag(sqrt(s(1:n)));
C=O(1:p,1:n);

% solve A from O by its shift property
A=O(1:end-p,:)\O(p+1:end,:);
% look for n that finds a stable model
E=sort(abs(eig(A)),'descend');
% semilogy(s,'-o');grid on,title(['n_{max} :' num2str(size(S,1))]) 
while ~isempty(find(E >=1))
   n=n-1; 
   O=U(:,1:n)*diag(sqrt(s(1:n)));
   C=O(1:p,1:n);   
   % solve A from O by its shift property
   A=O(1:end-p,:)\O(p+1:end,:);
   % look for n that finds a stable model
   E=sort(abs(eig(A)),'descend');
end
if n==0
    error('No stable A matrix could be found.')
end

% solve B,D from linear system
RHS=reshape(permute(reshape(U(:,n+1:end)'*L21*pinv(L11),[k*p-n,m^d,k]),[1,3,2]),[(k*p-n)*k,m^d]);
RHS(:,1)=0;
L=reshape(U(:,n+1:end)',[k*p-n,p,k]);
LHS=zeros(k*(k*p-n),p+n);
LHS(:,1:p)=reshape(permute(L,[1,3,2]),[k*(k*p-n),p]);
for i=2:k
   LHS((i-2)*(k*p-n)+1:(i-1)*(k*p-n),p+1:end)= reshape(L(:,:,i:end),[k*p-n,(k+1-i)*p])*O(1:(k-i+1)*p,:);
end
X=LHS\RHS;
D=X(1:p,:);
B=X(p+1:end,:);

end
