function order=armaorder(mvec,sig2,N,nu)
% Order estimation for a generic ARMA model using: AIC, AICc, GIC, and BIC.
% 
%  order=armaorder(mo,sig2,N,nu)
%       mvec    <- vector of model orders
%       sig2  <- vector mean square errors (that is, estimate of sigma^2) for model orders
%                given in mvec.
%       N     <- number of data points 
%       nu    <- GIC parameter  (usually nu \in [2,6]; default=4)
%       order -> the model order that minimizes the order estimation criterion.

% randy moses, 05 may 2005

if nargin<4
    nu=4;
end

% the general order estimation rule is -2 ln p_m + eta(m,N)*m
%
% where -2 ln p_m = N*ln(sigma^2_m) and where
%
%  for AIC:      eta(m,N) = 2 * m
%  for AIC_c:    eta(m,N) = 2 * (N)/(N-m-1) * m
%  for GIC:      eta(m,N) = nu * m 
%  for BIC(MDL): eta(m,N) = ln(N)* m

AIC  = N*log(sig2) + 2*mvec;
AICc = N*log(sig2) + (2*mvec * N)./(N-mvec-1);
GIC  = N*log(sig2) +  mvec * nu;
BIC  = N*log(sig2) +  mvec * log(N);

[~,b]=min(AIC);
order(1)=mvec(b);

[~,b]=min(AICc);
order(2)=mvec(b);

[~,b]=min(GIC);
order(3)=mvec(b);

[~,b]=min(BIC);
order(4)=mvec(b);