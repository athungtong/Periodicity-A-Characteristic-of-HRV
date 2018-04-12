function H=getentropy(p)
%H=getentropy compute entropy  -sum(p*log2(p))

if ~issparse(p)
    p=reshape(p,numel(p),1);
    p(p==0)=[];
else
    p=nonzeros(p);
end
H=-dot(p,log2(p));
    