function w=minnorm_r(r,n)
% The Root Min-Norm frequency estimator.
%
% w=minnorm(r,n);
%      r  ->  the true ACS
%      n  ->  the model order
%      w  <-  the frequency estimates

% Copyright 1996 by R. Moses

r=r(:);
m=length(r);

R=toeplitz(r);


% to use the forward-backward approach, uncomment the next line
% R=(R+fliplr(eye(m))*R.'*fliplr(eye(m)))/2;

% get the eigendecomposition of R; use svd because it sorts eigenvalues
[U,~]=svd(R);
S=U(:,1:n);
alpha = S(1,:)';
Sbar = S(2:m,:);

if norm(alpha) ~=1, 
   g = - Sbar * alpha / (1-alpha'*alpha);
else
   error('The min-norm solution does not exist');
end

% find the n roots of the a polynomial that are nearest the unit circle,
ra= conj(roots([1;g]));

% pick the n roots that are closest to the unit circle
[~,I]=sort(abs(abs(ra)-1));
w=angle(ra(I(1:n)));
