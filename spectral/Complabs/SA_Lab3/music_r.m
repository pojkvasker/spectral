function w=music_r(r,n)
% The Root MUSIC method for frequency estimation.
%
% w=music(r,n);
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
G=U(:,n+1:m);

C = G*G';

% find the coefficients of the polynomial in (4.5.16)
a = zeros(2*m-1,1);
for j=-(m-1):(m-1)
    a(j+m) = sum(diag(C,j));
end

% find the n roots of the a polynomial that are nearest and inside the unit circle,
ra=roots(a);
rb=ra(abs(ra)<1);

% pick the n roots that are closest to the unit circle
[~,I]=sort(abs(abs(rb)-1));
w=angle(rb(I(1:n)));
