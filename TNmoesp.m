function [A,BD,C]=TNmoesp(y,u,d,varargin)
% [A,BD,C]=TNmoesp(y,u,d) or [A,BD,C]=TNmoesp(y,u,d,n)
% --------------------------------------------------------
% Tensor network implementation of the MOESP algorithm for the identification of
% a specific polynomial state space model. 
%
% A 		= 	matrix, n x n state update matrix,
%
% BD 		=	tensor network structure, (p+n) x m^d matrix, m is the number of inputs,
%
% C 		=	matrix, p x n state output model matrix,
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
if k < 2
	error(['k is smaller than 2, please increase the sample size.'])
end
r=k*nchoosek(d+m-1,m-1)-(k-1);
N=r+k*p;

% construct block Hankel matrices
Y=zeros(k*p,N);
U=zeros(k*m,N);
for i=1:N
    U(:,i)=reshape(u(i:i+k-1,:)',[k*m,1]);
    Y(:,i)=reshape(y(i:i+k-1,:)',[k*p,1]);
end
U=reshape(U,[m,N*k]);
linvL11=rkh2tn(U,d);
linvL11.core{d}=reshape((linvL11.core{d}),[linvL11.n(end,1),m*k,N,1]);
linvL11.n(end,[2,3])=[m*k,N];

% Compute L21, L22 and (left) pseudo inverse of L11
if prod(linvL11.n(d,1:2)) >=linvL11.n(d,3)
    [U,S,V]=svd(reshape(linvL11.core{d},[prod(linvL11.n(d,1:2)),linvL11.n(d,3)]),'econ');
else
    error('Last core of U is underdetermined. Cannot do SVD. Input probably not PE.');
end
linvL11.core{d}=reshape(U,linvL11.n(d,:));
s=diag(S);
% inspect persistency of excitation of order k on input by checking the
% rank-gap
disp(['Check PE of input through rank gap of U matrix: ' num2str(s(r)/s(r+1))])
Sinv=[diag(1./s(1:r)),zeros(r,linvL11.n(d,3)-r)];
for i=1:d
   linvL11.core{i}=reshape(permute(reshape(linvL11.core{i},linvL11.n(i,:)),[1,3,2,4]),[linvL11.n(i,1),linvL11.n(i,3),linvL11.n(i,2),linvL11.n(i,4)]);
   linvL11.n(i,:)=[linvL11.n(i,1),linvL11.n(i,3),linvL11.n(i,2),linvL11.n(i,4)];
end
linvL11.core{d}=permute(reshape(Sinv*reshape(permute(linvL11.core{d},[2,3,1]),[linvL11.n(d,2),linvL11.n(d,3)*linvL11.n(d,1)]),[r,linvL11.n(d,3),linvL11.n(d,1)]),[3,1,2]);
linvL11.n(d,2)=r;
L=Y*V;
L21=L(:,1:r);
L22=L(:,r+1:end);

% determine system order n as numerical rank of L22
[U,S,~]=svd(L22,'econ');
s=diag(S)/S(1,1);
if ~isempty(varargin)
	if varargin{1} >= k
		disp(['Warning, given system order n should be strictly smaller than ' num2str(k) ' for this sample size.'])
		disp(['System order will be set to maximally ' num2str(k) '.']);
		n=k;
	else
	    n=min(varargin{1},k);
	end
else
    tol=numel(L22)*eps; % no output noise assumed   
    n=min(min(sum(s>tol),100),k);
end
O=U(:,1:n)*diag(sqrt(s(1:n)));
C=O(1:p,1:n);

% solve A from O by its shift property
A=O(1:end-p,:)\O(p+1:end,:);
% look for n that finds a stable model
E=sort(abs(eig(A)),'descend');
% semilogy(s,'-o');grid on,title(['n_{max} :' num2str(size(S,1))]) 
while ~isempty(find(E >=1,1))
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

% % solve B,D from linear system
L=reshape(U(:,n+1:end)',[k*p-n,p,k]);
LHS=zeros(k*(k*p-n),p+n);
LHS(:,1:p)=reshape(permute(L,[1,3,2]),[k*(k*p-n),p]);
for i=2:k
   LHS((i-2)*(k*p-n)+1:(i-1)*(k*p-n),p+1:end)= reshape(L(:,:,i:end),[k*p-n,(k+1-i)*p])*O(1:(k-i+1)*p,:);
end
% we overwrite the last core of linvL11 to form RHS 
BD=linvL11;
BD.core{d}=reshape(permute(reshape(U(:,n+1:end)'*L21*reshape(permute(linvL11.core{d},[2,3,1]),[r,linvL11.n(d,1)*linvL11.n(d,3)]),[k*p-n,m,k,linvL11.n(d,1)]),[4,1,3,2]),[linvL11.n(d,1),(k*p-n)*k,m,1]);
BD.n(d,2:3)=[k*(k*p-n),m];
% invert LHS and left-multiply with RHS
BD.core{d}=permute(reshape(pinv(LHS)*reshape(permute(BD.core{d},[2,3,1]),[(k*p-n)*k,m*BD.n(d,1)]),[n+p,m,BD.n(d,1)]),[3,1,2]);
BD.n(d,2)=n+p;

% set first column of BD to zero
% we contract with rank-2 mpo that represents 
% m^d x m^d Identity matrix with element (1,1) removed
I=eye(m);
for i=1:d
    if i==1
        temp=[I,zeros(m,m)];
        temp(1,m+1)=-1;
        temp=reshape(permute(BD.core{i},[4,1,2,3]),[BD.n(i,4)*prod(BD.n(i,1:2)),BD.n(i,3)])*temp;
        temp=reshape(temp,[BD.n(i,4),BD.n(i,1:3),2,1]);
        temp=permute(temp,[2,6,3,4,1,5]);
        BD.n(i,4)=BD.n(i,4)*2;
        BD.core{i}=reshape(temp,BD.n(i,:));
    elseif i==d
        temp=[I,zeros(m,m)];
        temp(1,m+1)=1;
        temp=reshape(permute(BD.core{i},[4,1,2,3]),[BD.n(i,4)*prod(BD.n(i,1:2)),BD.n(i,3)])*temp;
        temp=reshape(temp,[BD.n(i,4),BD.n(i,1:3),1,2]);
        temp=permute(temp,[2,6,3,4,1,5]);
        BD.n(i,1)=BD.n(i,1)*2;
        BD.core{i}=reshape(temp,BD.n(i,:));
    else
        temp=[I,zeros(m,3*m)];
        temp(1,3*m+1)=1;
        temp=reshape(permute(BD.core{i},[4,1,2,3]),[BD.n(i,4)*prod(BD.n(i,1:2)),BD.n(i,3)])*temp;
        temp=reshape(temp,[BD.n(i,4),BD.n(i,1:3),2,2]);
        temp=permute(temp,[2,6,3,4,1,5]);
        BD.n(i,[1,4])=BD.n(i,[1,4])*2;
        BD.core{i}=reshape(temp,BD.n(i,:));        
    end
end
end
