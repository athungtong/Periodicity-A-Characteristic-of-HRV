function [SampEn ApEn] = Entropy(norm,s,m,t,r,type)
% % [SampEn ApEn] = Entropy(norm,s,m,t,r) computes Sample Entropy and Approximate
% % Entropy when "maximum norm" is employed. SampEn is not defined when the correlation sum Cai and Cbi = 0 
% % (in this case SampEn may equal inf or NaN. This happen when, may be, r is too small.
% %
% % norm is the "maximum norm" of reconstructed embedding dimension (dim = m, delay t)
% % s is the time series.
% % r=tolorence, default =.2 (note: the time series should be normalized to
% % have unit varance
% % type = type of norm ('max' or 'Euc'). Default is 'max'
% % Reference: Richman, J. S., and Moorman, J. R. Physiological time-series 
% % analysis using approximate entropy and sample entropy. Am J Physiol Hreat Circ 
% % Physiol 278 (2000), 2039-2049
% %
% % This is the special case (for maximum norm only) of computing SampEn and ApEn. This function
% % speed up the computation time of the regular SampEn and ApEn function (about 10 times faster). 
% % The idea of speeding up the computation is the following. First, we
% know the way of how to compute the norm2 when m increase by one using
% knowledge of norm. By using idea, then we know that if any norm(i,j) such
% that r<norm(i,j), increasing dimension by one, will result in the new
% norm2(i',j') = max( norm(i,j), abs(something) ) which of course will
% again greater than r so it is not counted. So, the idea is we look only
% at all (i,j) such that r>norm(i,j) instead of all (i,j). This save lot of
% time in computing norm2. 
if nargin == 4
        r=.2; type='max';
end
if nargin == 5
    type='max';
end

Nb = (1+sqrt(1+8*length(norm)))/2; %(length(s)-(m-1)*t)
Na = Nb-t;%(length(s)-(m)*t)

index = find(r>=norm);
[I J] = get_embedded_index(index,Nb); 

%-------save this for ApEn
Bi = 1+histc([I J],1:Nb);
%-------------------------

% compute correlation sum of norm, which defind for i = 1:Nb-t,t=Nb-Na
I(J>Na)=[]; J(J>Na)=[];

if strcmp(type,'Euc')
    index=find(r>= norm((I-1)*Na-I.*(I+1)/2+(I-1)*t+J)+(s(I+m*t)-s(J+m*t)).^2);
else
    index=find(r>= abs(s(I+m*t)-s(J+m*t)));
end
Cai=2*length(index)/Na;
Cbi = 2*length(I)/Na;
SampEn = -log( Cai / Cbi);

% ApEn down here

Ai = 1+histc([I(index) J(index)],1:Na);
ApEn = 1/Nb*sum(log(Bi/Nb)) - 1/Na*sum(log(Ai/Na));
if isempty(index)
    ApEn=[];
end


% % % This is the spare code. May needed
% % 
% % % compute correlation sum of norm, which defind for i = 1:Nb-t,t=Nb-Na
% % index = find(r>=norm);
% % [I J] = get_embedded_index(index,Nb); 
% % 
% % %-------save this for ApEn
% % Bi = 1+histc([I J],1:Nb);
% % %-------------------------
% % % compute correlation sum of norm, which defind for i = 1:Nb-t,t=Nb-Na
% % I(J>Na)=[];
% % Cbi = 2*length(I)/Na;
% % 
% % [Cai Ai] = corrsum(norm2,r,0);
% % Ai=Ai+1;
% % SampEn = -log( Cai / Cbi);
% % ApEn = 1/Nb*sum(log(Bi/Nb)) - 1/Na*sum(log(Ai/Na));


