function [px nxbin atxbin] = get_margin_eq_prob(x,nxbin)
% [atxbin px] = get_margin_eq_prob(x,n) compute marginal entropy using
% marginal equinxbinuantization method. The marginal boxes are not defined
% equidistantly (in this case use get_margin_prob() ) but so that there is
% aproximately the same number of data points in each marginal bin. 
%
% x : input time series
% n : dimension of the data default = 1
% reference Milan Palus. see a book name: Information Theory and
% Statistical Learning page 197
%%
N = length(x);
if nargin ==1
    n=2;
    nxbin = round(N^(1/(n+1))); %number of bin
end


px = floor(N/nxbin)*ones(nxbin,1);
nrem = rem(N,nxbin);
ex=randperm(nxbin,nrem);
px(ex)=px(ex)+1;
platend=cumsum(px);
px = px/N;

platstart = [1 ; platend(1:end-1)+1];

[~, ix] = sort(x);

atxbin = zeros(N,1);
for i=1:length(platstart)
    atxbin(ix(platstart(i):platend(i))) = i;    
end

