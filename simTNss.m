function y=simTNss(A,BD,C,x0,u)
% y=simTNss(A,B,C,D,x0,u,d)
% ------------------------
% Simulates the polynomial state space system as specified by the A,C,
% D matrices, the BD tensor network, the initial state x0, the input
% vectors u and maximal total degree d.
%
% y         =   matrix, an L x p matrix where y(:,k) is the kth output,
%
% A 		= 	matrix, n x n state update matrix,
%
% BD 		=	tensor network structure, (p+n) x m^d matrix, m is the number of inputs,
%
% C 		=	matrix, p x n state output model matrix,
%
% x0		=	vector, n x 1 intial state vector,
%
% u 		=	matrix, L x m matrix where u(:,k) is the kth input,
%
% d 		=	scalar, maximal total degree of polynomial nonlinearity.
%
% Reference
% ---------
%
% Tensor network subspace identification of specific polynomial state space models
%
% 2017, Kim Batselier, Ching Yun Ko, Ngai Wong

p=size(C,1);
N=size(u,1);
y=zeros(N,size(C,1));
d=size(BD.n,1);
for i=1:N
    % contract [1,u(i,:)]' vector with BD TN
    yx=BD.core{1};
    for j=1:d-1
        yx=[1,u(i,:)]*reshape(yx,[BD.n(j,3),BD.n(j,4)]);
        yx=yx*reshape(BD.core{j+1},[BD.n(j+1,1),prod(BD.n(j+1,2:end))]);       
    end
    yx=reshape(yx,BD.n(end,2:3))*[1,u(i,:)]';    
    y(i,:)=(C*x0)+yx(1:p);
    x0=A*x0+yx(p+1:end);  
end

end
