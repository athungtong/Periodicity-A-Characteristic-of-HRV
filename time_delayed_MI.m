function [MI optT] = time_delayed_MI(x,tao)
% MI = time_delayed_MI(x,tao) computes mutual information between x(t) and
% x(t+tao),tao must be column vector has value between [0 length(x)-1]
% optT is the time index corresponding to the
% first minimum MI used to find optimum Time delay in embbedding
% dimension
% Reference: A. M. Fraser, and H. L. Swinney, 
% Independent coordinates for strange attractors from mutual information, 
% Phys. Rev. A 33 (1986) 1134

if nargin == 1
    tao=0:round(length(x)/4);
end

T=length(tao);
MI=zeros(T,1);
for t=0:T-1
    MI(t+1)=mutual_info(x(1:end-tao(t+1)),x(1+tao(t+1):end),1);
end

switch nargout
    case{0}
        plot(tao,MI); ylabel('mutual information');xlabel('time delayed');
    case{2}
        [Tmax Tmin]=getlocalmaxmin(MI);
        optT=tao(find(Tmin==1, 1 )); %find first minimum
end
