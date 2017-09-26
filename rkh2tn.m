function TN=rkh2tn(U,varargin)
% TN=rkh2tn(U) or TN=rkh2tn(U,d)
% ------------------------------
% When given a single cell argument U, then the Tensor Network corresponding
% with the matrix formed from the Khatri-Rao product of U{1},U{2},...,U{d}
% is formed. If U is a matrix then the Tensor Network of the d-times repeated
% Khatri-Rao  product of U is formed.
%
% TN 	=   Tensor Network structure,
%
% U 	=   cell/matrix,
%
% d 	=	scalar, argument only required for matrix U to indicated how many times
%           its Khatri-Rao product needs to be taken.
%
% Reference
% ---------
%
% 09/2017, Kim Batselier

if iscell(U)
    d=length(U);    
else
    d=varargin{1};
    temp=U;
    U=cell(1,d);
    for i=1:d
        U{i}=temp;
    end
end

n=zeros(d,1);
for i=1:d
    % second dimension has to be the same over all factor matrices    
    [n(i),N]=size(U{i});    
end
TN.core=cell(1,d);
TN.n=ones(d,4);

% initialize first core
TN.n(1,:)=[1,n(d),N,1];
TN.core{1}=reshape(U{1},TN.n(1,:));

for i=1:d-1
    temp=reshape(TN.core{i},[prod(TN.n(i,1:2)),prod(TN.n(i,3:4))]);
    temp=khatr(temp,U{i+1});
    temp=reshape(temp,[TN.n(i,1)*n(i),n(i+1)*N]);
    [U1,S1,V1]=svd(temp,'econ');
    s=diag(S1);
    tol=eps(s(1))*max([N,n(i),n(i+1)]);
    r=sum(s>tol);
    TN.n(i,2:end)=[n(i),1,r];
    TN.core{i}=reshape(U1(:,1:r),TN.n(i,:));
    TN.n(i+1,1:end-1)=[r,n(i+1),N];
    TN.core{i+1}=reshape(S1(1:r,1:r)*V1(:,1:r)',TN.n(i+1,:));
end
end
