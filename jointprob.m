function [H p]=jointprob(x,y,z,N,atxbin,nxbin,atybin,nybin,atzbin,nzbin)
% [H p]=jointprob computes joint probability of 2 or 3 variables
% usage: 
% H = jointprob(x,y) need only Entropy (-dot(p,log2(p)), this is much
% faster than getting p 

% [H p] = jointprob(x,y)        :2 variables also need p
% p = jointprob(x,y,z)      :3 variables
% p = jointprob([],[],[],N,atxbin,nxbin,atybin,nybin)    :2 variables
% know atxbin, N = length of x
% p = jointprob([],[],[],N,atxbin,nxbin,atybin,nybin,atzbin,nzbin) :3 variables know atxbin
% NOTE: p nxbin*nybin*nzbin x 1 sparse vector

numvar=0;
if nargin == 2 || nargin == 3
    N=size(x,1);
    nxbin=16;%round(sqrt(N));
    nybin=nxbin; 

    [~, ~, atxbin]= get_margin_eq_prob(x,nxbin);
    [~, ~, atybin]= get_margin_eq_prob(y,nybin);
    numvar=2;
end

if nargin==3
   nzbin=nxbin;
   [~, ~, atzbin]= get_margin_eq_prob(z,nzbin);
   numvar=3;
end

% use sparse to reduce memory
if nargin == 8 || numvar==2     
    H=0;
    p = sparse(nxbin*nybin,1);
    numr=[0 nxbin:nxbin:nxbin*nybin];
    for i=1:nxbin
        indx=find(atxbin==i);
        if isempty(indx),continue; end
        atybin_i = atybin(indx);
        for j=1:nybin
            nij=length(find(atybin_i==j));
            if nij == 0, continue; end
            H=H-nij/N*log2(nij/N);
            if nargout==2
                vecind=numr(j)+i;
                p(vecind)=nij;
            end
        end
    end
    
elseif nargin==10 || numvar==3
    H=0;
    p = sparse(nxbin*nybin*nzbin,1);
    numr=[0 nxbin:nxbin:nxbin*nybin];
    nump = [0 nxbin*nybin:nxbin*nybin:nxbin*nybin*nzbin];
    for i=1:nxbin
        indx=find(atxbin==i);
        if isempty(indx),continue; end
        atybin2=atybin(indx);        atzbin2=atzbin(indx);
        for j=1:nybin
            indy=find(atybin2==j);
            if isempty(indy),continue;end
            atzbin3=atzbin2(indy);
            for k=1:nzbin
                nijk=length(find(atzbin3==k));
                if nijk == 0,continue;end
                H=H-nijk/N*log2(nijk/N);
                if nargout==2
                    vecind = nump(k)+numr(j)+i;
                    p(vecind)=nijk;
                end
            end
        end
    end
    
else
    warning('inappropriate number of input argument');
    p=[];
end

p=p/N; %p is sparse matrix