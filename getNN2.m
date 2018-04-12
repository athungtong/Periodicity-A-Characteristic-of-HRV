function [rnn knn] = getNN2(norm,k,P)
%%
% rnn = getNN2(norm,k, P) returns nearest neighbor of each embedded point i for
% i = 1,2,...,Ne-k with a constrain that |i-j|>P. The constrain avoid taking
% the adjacent points to be the nearest neighber. In
% this case, we consider only points on other loop. T can be the mean or
% median period of the time series. norm is the Euclidian or maximum norm 
% Uncomment some code below to get the index together with the min value

% Input Arguments
% norm is the Euclidean distance between point (i,j), i<j in the embedding space with
%      dimension m, and time delay t (norm = getnorm(x,m,t)). 
% k. The number of point you like to excluded. This is useful if you like
% to compare rnn at dimension m and rnn at dimension m+1. If you set k = t,
% then dimension of rnn(m) would be t smaller so that it is equal to
% dimension of rnn(m+1) and then you can compare.
% P is the period (defind as number of point) of the time series
%
% Output Arguments
% rnn: distance from points i, all i, to its nearest neighbor
% knn: index "in norm vector domain" indicating which point is cloest to
% which point. The index in vector is the following
%
% 12 13 14 15
%    23 24 25
%       34 35
%          45
% If point 1 is closest to point 5, then knn(1) = 4. If point 2 is closest to point 1 then
% knn(2)=1. and so on. You can convert knn to the index in the vector space
% using the command
%
% [I J] = get_embedded_index(knn,Ne);
% For example, if knn = 4 then the output I=1 and J=5 which is the pair
% that we found above
%
switch nargin
    case{1}
        k=0;P=0;
    case{2}
        P=0;
end

Ne = (1+sqrt(1+8*length(norm)))/2; %number of reconstructed points
J=(1:Ne)*Ne-(1:Ne).*((1:Ne)+1)/2;

rnn=zeros(Ne-k,1);

if nargout==1
    rnn(1) = min(norm((1+P):Ne-1-k));
    for i = 2:Ne-k
        rnn(i) = min(norm([J(1:i-1-P)-(Ne-i) J(i-1)+1+P:J(i)-k]));
    end
else
    knn=zeros(1,Ne-k);    
    I=(1+P):Ne-1-k;
    [rnn(1) inn] = min(norm(I));
    knn(1) = I(inn);
    for i = 2:Ne-k
        I=J(1:i-1-P)-(Ne-i);
        [rnn(i) inn] = min(norm(I));
        knn(i) = I(inn);
    end
end