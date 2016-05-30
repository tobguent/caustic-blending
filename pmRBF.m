function [V,f] = pmRBF(xc,x,f)
% evaluate radial basis functions f in d-dimensional domain
% - xc centers n x d
% - x evalution points m x d
% - f function handle, default is thin plate splines f(r2)=r2*log(sqrt(r2))
%     with squared radius r2=||xc-x||^2
% Result
% - V matrix of values m x n

if nargin<3 || isempty(f),
    f=@(r2)(r2.*log(sqrt(r2)));
end

m=size(x,1);
n=size(xc,1);
d=size(xc,2);

assert(d==size(x,2));

x = reshape(x',d,1,[]);

r2 = repmat(xc',[1,1,m]) - repmat(x,[1,n,1]);
r2 = sum(r2.^2,1);

phi = f(r2);
phi(~isfinite(phi))=0; % e.g., log(0) at center

V = reshape(reshape(phi,1,[]),n,[])';
