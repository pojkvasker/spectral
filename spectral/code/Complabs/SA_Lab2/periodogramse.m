function phi=periodogramse(y,v,L)
% The windowed periodogram spectral estimator for estimating the power spectral density (PSD).
%
% phi=periodogramse(y,v,L)
%      y   -> the data vector
%      v   -> the window vector of length M
%      L   -> the number of PSD samples
%      phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

if ~isvector(y)
    error('The input must be a data vector, not a matrix.')
end

% check the length of the window vector
M=length(v);
N=length(y);
if (M>N)
   error('The length of the window is larger than the length of the data vector');
elseif (M<N)
   warning('The length of the window is smaller than the length of the data vector; the data vector will be truncated to the window length')
end

y=y(:);         % columlize the data matrix 
% generate the spectral estimate

phi=(abs(fft((y(1:M).*v(:)),L)).^2)/M;
