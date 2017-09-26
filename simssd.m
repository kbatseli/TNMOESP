function y=simssd(A,B,C,D,x0,u,d)
% y=simssd(A,B,C,D,x0,u,d)
% ------------------------
% Simulates the polynomial state space system as specified by the A,B,C,
% D matrices, the initial state x0, the input vectors u and maximal total degree d.
%
% y         =   matrix, an L x p matrix where y(:,k) is the kth output,
%
% A 		= 	matrix, n x n state update matrix,
%
% B 		=	matrix, n x m^d input-to-state matrix, m is the number of inputs,
%
% C 		=	matrix, p x n state output model matrix,
%
% D 		= 	matrix, p x m^d input-to-output matrix,
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

N=size(u,1);
y=zeros(N,size(C,1));
for i=1:N
    U=mkron([1,u(i,:)]',d);
    y(i,:)=(C*x0)+D*U;
    x0=A*x0+B*U;  
end

end
