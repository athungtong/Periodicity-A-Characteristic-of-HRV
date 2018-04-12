function maxnorm = getmaxnorm(s,m,t)
% maxnorm = getmaxnorm(s,m,t) find infinite norm (maximum norm) of embedded s with
% dimension m (default=2) and time delay t (default=1)
switch nargin
    case{1}
        m=2;t=1;
    case{2}
        t=1;
end

if size(s,1)==1
    s=s';
end

Ne = length(s)-(m-1)*t;
maxnorm = zeros((Ne-1)*Ne/2,1);
% Special when m==1
if m == 1
    for i=1:Ne-1
        n=(i-1)*Ne-i*(i+1)/2;
        maxnorm(n+(i+1:Ne)) = abs(s(i)-s(i+1:Ne));
    end
    return;
end

I = zeros(m,Ne-1);
I(1,:)=(2:Ne);
for l=2:m
    I(l,:) = I(l-1,:)+t;
end

for i = 1:Ne-1
    n=(i-1)*Ne-i*(i+1)/2;
    maxnorm(n+(i+1:Ne))=max(abs(s(i:t:i+(m-1)*t)*ones(1,Ne-i)-s(I(:,i:Ne-1))));    
end    

