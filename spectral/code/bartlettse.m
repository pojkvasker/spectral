function phi=bartlettse(y,M,L)
% The Bartlett method for estimating the power spectral density (PSD).
%
% phi=bartlettse(y,M,L)
%      y -> the data vector
%      M -> the length of subsequences of y
%      L -> the number of PSD samples
%      phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

if ~isvector(y)
    error('The input must be a data vector, not a matrix.')
end

% check the length M
N=length(y);
if (M>N)
   error('M is greater than the data length.');
end

phi=welchse(y,ones(M,1),M,L);   % bartlett is a special case of welch.
