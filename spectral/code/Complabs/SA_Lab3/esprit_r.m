function w=esprit_r(r,n)
% The ESPRIT method for frequency estimation.
%
% w=esprit(r,n);
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
[U,D,V]=svd(R);
S=U(:,1:n);

phi = S(1:m-1,:)\S(2:m,:);

w=-angle(eig(phi));
