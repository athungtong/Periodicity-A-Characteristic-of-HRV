
function [I J] = get_embedded_index(N,Ne)
% [I J] = get_embedded_index(N,Ne) maps index in norm vector (computed by
% getnorm()) to the corresponding point (i,j),i<j in the embedding space. 
% N is the scalar or vector if norm index
% Ne is the size of embedding space (NexNe)
% [I J] is the scalar (if N is scalar) or row vector (size = 1xlength(N)) 

if size(N,1)>1
    N=N';
end

I = ceil((Ne-.5)-sqrt( (Ne-1/2)^2-2*N));
if nargout == 2
    J = Ne- I*Ne + I.*(I+1)/2 + N;
end
