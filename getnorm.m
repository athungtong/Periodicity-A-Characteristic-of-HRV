function norm = getnorm(s,m,t)
% norm = getnorm(s,m,t) compute Euclidean^2 norm of s based on embedded
% dimension m (default=2) and time delay t (default =1). Note that if you
% need to compute the Euclidean norm, you need to do the sqrt(norm).
% Remark, this function may suffer from numerical error. To avoid this use
% the very basic one below.
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
norm = zeros((Ne-1)*Ne/2,1);

%Special case for m = 1;
if m==1
    for i=1:Ne-1
        n=(i-1)*Ne-i*(i+1)/2;
        norm(n+(i+1:Ne)) = (s(i)-s(i+1:Ne)).^2;
    end
    return;
end

% This is special case for m=2.
if m==2
    for i = 1:Ne-1
        n=(i-1)*Ne-i*(i+1)/2;
        norm(n+(i+1:Ne)) = (s(i)-s(i+1:Ne)).^2 + (s(i+t)-s((i+1:Ne)+t)).^2;
    end
    return;
end

% This is for m not equal to 2
I = zeros(m,Ne-1);
I(1,:)=(2:Ne);
for l=2:m
    I(l,:) = I(l-1,:)+t;
end

for i = 1:t
    n=(i-1)*Ne-i*(i+1)/2;
    norm(n+(i+1:Ne))=sum( ( s(i:t:i+(m-1)*t)*ones(1,Ne-i) - s(I(:,i:Ne-1)) ).^2 );
end

for i = t+1:Ne-1
%This is to help getting index easier and faster
    n=(i-1)*Ne-i*(i+1)/2;
    c=(i-t-1)*Ne-(i-t-1)*( i-t)/2-i;
    norm(n+(i+1:Ne)) = ...
    norm(c+(i+1:Ne)) + (s(i+(m-1)*t)-s((i+1:Ne)+(m-1)*t)).^2-(s(i-t)-s((i+1:Ne)-t)).^2;
end


% % This is the simpler version, save normR as NexNe matrix
% % It is easier to call index but it is way slower than above
% Ne = length(s)-(m-1)*t;
% normR = sparse(Ne,Ne);
% for i = 1:t
%     si = s(i:t:i+(m-1)*t);
%     for j = i+1:Ne
%             normR(i,j) = sum((si-s(j:t:j+(m-1)*t)).^2);
%     end
% end
% 
% for i = t+1:Ne-1
%     si = s(i+(m-1)*t); sii = s(i-t);
%     for j = i+1:Ne
%         normR(i,j) = normR(i-t,j-t)+(si-s(j+(m-1)*t))^2 - (sii-s(j-t))^2;
%     end
% end

% % Also, this is the very basic version. This version avoid recursive method.  
% Ne = length(s)-(m-1)*t;
% norm = zeros((Ne-1)*Ne/2,1);
% I = zeros(m,Ne-1);
% I(1,:)=(2:Ne);
% for l=2:m
%     I(l,:) = I(l-1,:)+t;
% end
% 
% for i = 1:Ne-1
%     n=(i-1)*Ne-i*(i+1)/2;
%     norm(n+(i+1:Ne))=sum((s(i:t:i+(m-1)*t)*ones(1,Ne-i)-s(I(:,i:Ne-1))).^2);
% end


