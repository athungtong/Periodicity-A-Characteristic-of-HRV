function M = degreeofmono(s)
% compute degree of mono of s. If s is a vector, the degree of mono will be
% computed based on each eliment of s. If s is a matrix, the degree of mono
% will be compute based on each column of s. In other word, it will compare
% distribution of each column of s, let call s1,s2,s3... the comparison
% will based on rank sum test.

if min(size(s)) == 1 %if s is a vector
    r = length(s);
    M = 0;
    for i = 1:r-1
        for j = i+1:r
            M = M+sign(s(j)-s(i));
        end
    end
    M = 2*M/r/(r-1);
else % if s is a matrix, we compare distribution instead of a value
    r=size(s,2); M = 0;
    for i = 1:r-1
        for j=i+1:r
            [p,h,stats]=ranksum(s(:,j),s(:,i),'alpha',0.05);
            M = M+h*sign(stats.zval);
        end
    end
    M = 2*M/r/(r-1);
end