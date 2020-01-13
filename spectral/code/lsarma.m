function [a,b,sig2]=lsarma(y,n,m,K)
% The two-stage Least-Squares ARMA method given in section (3.7.2)
%
% [a,b,sig2]=lsarma(y,n,m,K)
%      y     -> the data vector
%      n     -> AR model order
%      m     -> MA model order
%      K     -> the order of the truncated AR model 
%      a     <- the AR coefficient vector
%      b     <- the MA coefficient vector
%      sig2  <- the noise variance

% Copyright 1996 by R. Moses

if ~isvector(y)
    error('The input must be a data vector, not a matrix.')
end

y=y(:);
N=length(y);             % data length
L=K+m;                   % This choice is one out of many possible

% N-L should be >= n+m (i.e. K<=N-n-2m)
if (N-L) < n+m || K>=N/2-1,
  error('K is too large');
end

% L>max(n,m) needed in general, here K>=n-m is sufficient as L>m by choice
if K<n-m
  error('K is too small');
end

% estimate alpha coefficients 
alpha=lsar(y,K);
  
% estimate the noise sequence e(t)
e=filter(alpha,1,y);       

% construct the z vector and Z matrix in equations (3.7.12) and (3.7.13)
z=[y(L+1:N)];
Z=[toeplitz(y(L:N-1),y(L:-1:L-n+1).'),-1*toeplitz(e(L:N-1),e(L:-1:K+1).')];

% estimate theta (a,b)
%theta=-Z\z;
theta=-pinv(Z)*z;
a=[1;theta(1:n)];
b=[1;theta(n+1:m+n)];

% estimate the noise covariance
sig2=norm(Z*theta+z)^2/(N-L);
