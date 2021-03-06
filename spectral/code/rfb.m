function phi=rfb(y,K,L)
% The statically stable Refined Filter Bank spectral estimator estimating 
% the power spectral density (PSD).
%
% phi=rfb(y,K,L)
%    y   <- the data vector (length N)
%    K   <- the ratio of the baseband filter bandwidth to 1/N (K=N*Beta)
%    L   <- the number of estimated spectral samples
%    phi -> the estimated spectrum

% Copyright 1996 by R. Moses

if ~isvector(y)
    error('The input must be a data vector, not a matrix.')
end

y=y(:);
N=length(y);       % data length

h=slepian(N,K,K);  % get the first K Slepian filters

phi=abs(fft(flipud(h).*(y*ones(1,K)),L)).^2;
if (K > 1)
   phi=mean(phi')';
end
