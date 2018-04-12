function norm2 = get_norm_mplusone(s,m,t,norm,type)
% % norm2 = get_norm_mplusone(s,m,t,norm,type) you can use this fnct
% %to find norm of time series when m increase by 1 but t and s is the
% %same. input type = 'Euc' if norm is euclidean norm (default is type='max')

% we can use information from the original m dimension. For norm(i),
% %we just need to add the different of the last point as the following. Note
% %that for dimension m, Ne = length(s)-m*t+t, for m+1,
% %Ne=length(s)-m*t, so the new one has t length less than the previous one
% %so we need to be careful about the index. Following is the code to find
% %norm of dimension m+1
if nargin == 4
    type='max';
end

if size(s,1)==1
    s=s';
end

Ne = length(s)-(m-1)*t; 
norm2 = zeros((Ne-1)*Ne/2,1);

%Compute Euclidean norm here.
if strcmp(type,'Euc')
    for i=1:Ne-1
        n=(i-1)*Ne-i*(i+1)/2;
        c=n+(i-1)*t;
        norm2(n+(i+1:Ne)) = norm(c+(i+1:Ne)) + (s(i+(m-1)*t)-s((i+1:Ne)+(m-1)*t)).^2;    
    end
else   %Compute maximum norm
    for i=1:Ne-1
        n=(i-1)*Ne-i*(i+1)/2;
        c=n+(i-1)*t;
           norm2(n+(i+1:Ne))= max( norm(c+(i+1:Ne)) , abs(s(i+(m-1)*t)-s((i+1:Ne)+(m-1)*t))); %for m=2
    end
end
end
