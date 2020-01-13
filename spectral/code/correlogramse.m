function phi=correlogramse(y,L)
% The correlogram spectral estimator for estimating the power spectral density (PSD).
%
% phi=correlogramse(y,L)
%    y   -> the data vector
%    L   -> the number of psd samples
%    phi <- spectral estimates at L frequencies w=0, 2*pi/L, ..., 2*pi(L-1)/L

% Copyright 1996 by R. Moses

if ~isvector(y)
    error('The input must be a data vector, not a matrix.')
end

y=y(:);         % columlize the data matrix
phi=periodogramse(y,ones(size(y)),L);  % the correlogram SE is the same as 
                                       % the unwindowed periodogram SE
