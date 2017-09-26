function y=khatr(varargin)
% y=khatr(A,B) or y=khatr(A,B,C)
% ----------------------------------
% Computes the column-wise right-Kronecker product of matrices A,B,C. Note
% that this computes the right-Kronecker product such that the order of the
% indices is maintained.
%
% y,A,B,C   = matrix.
%
% Reference
% ---------
%
% Batselier Kim, adapted from Zhongming's dotkron

if nargin==2
    L=varargin{1};     R=varargin{2};
    [r1,c1]=size(L);   [r2,c2]=size(R);
    if c1 ~= c2
        error('Matrices should have equal columns!');
    else
        y=repmat(L,r2,1).*kron(R,ones(r1,1));  % faster than using FOR
    end
elseif nargin==3
    L=varargin{1};     M=varargin{2};     R=varargin{3};
    c1=size(L,2);      c2=size(M,2);      c3=size(R,2);
    if c1 ~= c2 || c2 ~=c3  ||  c1~=c3
        error('Matrices should have equal columns!');
    else
        y=khatr(L, khatr(M, R));
    end
else
    error('Please input 2 or 3 matrices!');
end
    


end