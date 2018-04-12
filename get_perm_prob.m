% function [px xbin atxbin]=get_perm_prob(x,ord,permlist) compute permutation probability of
% time series x with order (ord). permlist is perms(1:ord) which is all
% possible permutation set of number 1,2,3,...,ord. if permlist is not
% given, it will be computed inside this function. 
% 
% Reference: Permutation Entropy: A Natural Complexity Measure for Time Series
% Bandt, Christoph
% Pompe, Bernd
% VOLUME 88, NUMBER 17 PHYSICALREVIEW LETTERS 29 APRIL 2002

function [px nxbin atxbin]=get_perm_prob(x,ord,xbin)

if nargin<3
    xbin=perms(1:ord);
end
nxbin = size(xbin,1);

%need column vector
if size(x,1)>size(x,2),x=x';end

N=length(x);
% clear c
c=zeros(nxbin,1);
atxbin = zeros(N-ord+1,1);
[~,permlist]=sort(x(1:ord)); 
binlist = find(ismember(xbin,permlist,'rows'));
atxbin(1)=binlist;
c(binlist)=1;

for j=2:N-ord+1
    [~,iv]=sort(x(j:j+ord-1)); 
    newlist=1;
    for i=1:size(permlist,1)
        if permlist(i,:)==iv
            atxbin(j)=binlist(i);
            c(binlist(i)) = c(binlist(i))+1;
            newlist=0;
            break
        end
    end
    if newlist        
        permlist=[permlist;iv];
        temp = find(ismember(xbin,iv,'rows'));
        binlist=[binlist;temp];
        atxbin(j) = temp;
        c(temp)=1;
    end   
end

px=c/(N-ord+1); 

