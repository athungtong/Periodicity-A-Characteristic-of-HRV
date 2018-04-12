function [gamma t]= acf(x,n,c)
%%
% gamma=acf(x,n,c) plot sample autocorrelation of x at delay 1 to n. If n is
% not defined, the default is half of the length of x. plot color is c,
% default is blue

x=center(x);
% x=x-mean(x);
if nargin == 1
    n=round(length(x)/2); c='b';
end
if nargin==2
    c='b';
end
N = length(x);
mx = mean(x);
gamma = zeros(2*n+1,1);
gamma(n+1) = sum( (x-mx).*conj(x-mx))/(N) ;
for t = 1:n
    gamma(n+1+t) = sum( (x(1:N-t)-mx).*conj(x(1+t:N)-mx))/(N-t) ;
    gamma(n+1-t) = gamma(n+1+t);
end
t = -n:n;
if nargout == 0
    plot( t, gamma ,c);
end

end