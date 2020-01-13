function [beta,relMSE,yHat]=lsa(y,w)
% Computes the complex amplitudes of sinusoidal components given the 
% frequencies by a least-squares fit to the data.
% Model:
%   y(t)=sum_{k=1}^n  beta_k e^{i w_k t} + e(t),   t= 0,...,N-1 
%
% [beta,relMSE,yHat]=lsa(y,w)
%       y       -> Nx1 data vector
%       w       -> sinusoidal frequencies
%       beta    <- the complex amplitude 
%       rel_MSE <- the relative error between the data and its reconstruction
%       yHat    <- reconstructed signal using the frequencies and amplitudes

y=y(:);
w=w(:);

%Vandermonde matrix of frequencies
A=exp(1j*(0:length(y)-1)'*w');
%Compute amplitudes using LS
beta=A\y;
%Reconstruct the signal
yHat=A*beta;

%Ensure that the result is real valued if y is real valued (numerical errors can cause complex result)
if isreal(y) 
    yHat=real(yHat);
end

%Compute the relative MSE
relMSE=1-yHat'*yHat/(y'*y);

