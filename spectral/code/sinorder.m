function order=sinorder(mvec,sig2,N,nu)
% Order estimation for sinusodial models using: AIC, AICc, GIC, and BIC.
% 
% order=sinorder(mvec,mse,N,nu)
%       mvec  <- vector of number of sinusoids (or complex exponentials for complex valued data)
%       sig2   <- vector mean square errors (that is, estimate of sigma^2) for model orders
%                given in mvec.
%       N     <- number of real-valued data points 
%       nu    <- GIC parameter  (usually nu \in [2,6]; default=4)
%       order -> the model orders that minimizes the AIC, AICc, GIC, and BIC criterions.

% randy moses, 09 sep 2003

if nargin<4, 
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
%
% for the first three methods, m=3*n_c+1, where n_c=#sinusoids
% for BIC, m=5*n_c+1, where n_c=#sinusoids

AIC  = N*log(sig2) + 2*(3*mvec+1);
AICc = N*log(sig2) + (2*(3*mvec+1) * N)./(N-(3*mvec+1)-1);
GIC  = N*log(sig2) +  (3*mvec+1) * nu;
BIC  = N*log(sig2) +  (5*mvec+1) * log(N);

[~,b]=min(AIC);
order(1)=mvec(b);

[~,b]=min(AICc);
order(2)=mvec(b);

[~,b]=min(GIC);
order(3)=mvec(b);

[~,b]=min(BIC);
order(4)=mvec(b);
return;
