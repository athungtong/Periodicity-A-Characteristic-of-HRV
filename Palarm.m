%   Palarm.m
function [Out,Kout]=Palarm(X,pf,m)
%%

%   Palarm.m
%   usage: [Out,Kout]=Palarm(X,pf)
%   Detects change-points in the mean of a random sequence
%
%   Anatoly Zlotnik and Alexandra Piryatinska, 2004
%   Revised 2005 by Anatoly Zlotnik
%
%   Input:
%
%   X       time-series data
%   pf      false alarm probability 
%
%   Notes:  Algorithm details:
%   Parameter   Default Description
%   delta1      1       weighting parameter for phase I
%   delta2      0       weighting parameter for phase II
%   pf          0.1     false detection probability
%   m           3       minimum interval size
%   M           1000    number of terms in KS series
%   epsilon     0.02    minimum relative distance of change point from 
%                       endpoints of local interval
%   
%   Output:
%
%   Out                 data with each interval replaced by local mean
%   Kout                vector containing estimated global change points
%  
%   Anatoly Zlotnik, June 19, 2006
%
%   [1] Zlotnik, A., "Algorithm Development for Estimation and Modeling 
%       Problems In Human EEG Analysis". M.S. Thesis, Case Western Reserve 
%       University; 2006
%

%Estimation of change points and mean on each stable interval 
warning off MATLAB:fzero:UndeterminedSyntax
if(size(X,2)>1), X=X'; end
if(size(X,2)>1), disp('error in input!'); return; end
%default parameters
% if(nargin<7) epsilon=0.02; end
% if(nargin<6) M=1000; end; 
% if(nargin<5) m=100; end
% if(nargin<4) delta2=0; end; 
% if(nargin<3) delta1=1; end
% if(nargin<2) pf=0.1; end;

if(nargin<2)
    pf=0.1;m=3;
end
if nargin<3
    m=3;
end
delta1=1; delta2=1;  delta3=1; M=1000; epsilon=0.02;

%scale input
Xin=X; Xbar=mean(X); X=(X-Xbar)/std(X);
N=length(X);

% calculate forward and backward cumulate sum of X. These variable will be
% used in computing family statistic ystat. They will extremely speed up
% the computation.

CX1=cumsum(X);
CX2=flipud(cumsum(flipud(X(2:end))));
CXX=cumsum(X.^2);

%Calculate Initial change point threshold
thresh=abs(KSinverse(pf,M));

%Calculate initial change points (Phase I)
P0=1; Kin0=[]; t1=1;t2=N;
[Kone]=PalarmF(CX1,CX2,CXX,t1,t2,Kin0,P0,delta1,thresh,m,epsilon);
if(isempty(Kone)), Kout=Kone; Out=mean(Xin)*ones(N,1); return; end
Kone=sort(Kone);
%Calculate new threshold
thresh=abs(KSinverse(pf/10,M));

% plot(Kone,'g+');hold on;

%Eliminate doubtful change points (Phase II)
[Kout]=diagn(X,CX1,CX2,CXX,Kone,delta2,m,thresh);
% plot(Kout,'g+');hold on;

%Estimate final locations of change points (Phase III)
Kout=finalal(CX1,CX2,CXX,Kout,m,delta3);

% [Kout]=diagn2(CX1,CX2,CXX,Kout,delta2,m,thresh,epsilon);
% [Kout]=diagn(X,CX1,CX2,CXX,Kout,delta2,m,thresh);
% 
% Kout=finalal(CX1,CX2,CXX,Kout,m,delta3);
% Kout=finalal(CX1,CX2,CXX,Kout,m,delta3);

% plot(Kout,'g+');hold on;
Out=getmeanK(Xin,Kout);

if(size(Out,2)>1), Out=Out'; end

function [Kout]=PalarmF(CX1,CX2,CXX,t1,t2,Kin,P,delta,thresh,m,epsilon)
%[Kout]=PalarmF(X,Kin,P,delta,thresh,m,epsilon)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   X       local data
%   Kin     vector containing estimated global change points
%   P       global index of first point of local data
%   delta   detect parameter
%   thresh  minimum statistic for change point characterization
%   m       minimum interval size
%   epsilon minimum relative distance of change point from endpoints
%
%   Kout    updated vector of estimated global change points
%
%Recursive algorithm for detection of change points by successive interval
%bisection

[stat,k,Sigma]=ystat(CX1,CX2,CXX,t1,t2,delta,m);
N=t2-t1+1;%length(Xin);
% Sigma=std([(Xin(1:k)-mean(Xin(1:k))) (Xin(k+1:N)-mean(Xin(k+1:N)))])/sqrt(N);

Level=thresh*Sigma;
D=max([round(epsilon*N) 1]);
if(stat(k)<=Level), Kout=Kin; return;
else
%     X1=Xin(1:k-D);
    P1=P; t11=t1; t12=t1+k-D-1;
    
%     X2=Xin(k+D:N);
    P2=P+k+D-1; t21=t1+k+D-1;   t22=t1+N-1;
    
    if k-D <= m %if(length(X1)<=m)
        Ktemp1=Kin; 
    else
        Ktemp1=PalarmF(CX1,CX2,CXX,t11,t12,Kin,P1,delta,thresh,m,epsilon); 
    end
    
    if N-(k+D)+1 <= m %if(length(X2)<=m)
        Ktemp2=Ktemp1; 
    else
        Ktemp2=PalarmF(CX1,CX2,CXX,t21,t22,Ktemp1,P2,delta,thresh,m,epsilon); 
    end
    
    if k-D > m && N-(k+D)+1 > m %if(length(X1)>m && length(X2)>m)
        Kout=[Ktemp2; P+k-1]; 
    else
        Kout=Ktemp2; 
    end
end

function [ystat1,k,Sigma] = ystat(CX1,CX2,CXX,t1,t2,delta,m)
%[ystat1,k,Sigma]=ystat(X,delta,m)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   X       Input data
%   delta   detect parameter
%   m       minimum interval size
%
%   ystat   test statistic
%   k       estimated change point
%
%calculate statistics for change point and maximize
%to detect change point k

t=(1:(t2-t1))'; N=t2-t1+1;
tflip = ((t2-t1):-1:1)';

if(N<m), ystat1=zeros(N,1); k=1; return; end

% M1=cumsum(X(1:N-1))./[1:N-1]'; 
% M2=flipdim(cumsum(flipdim(X(2:N),1))./[1:N-1]',1);

if t1==1
    M1=CX1(t1:(t2-1));
    MXX=CXX(t2);
else
    M1=(CX1(t1:(t2-1))-CX1(t1-1));
    MXX=CXX(t2)-CXX(t1-1);
end
if t2==length(CX1)
    M2=CX2(t1:(t2-1));
else
    M2=(CX2(t1:(t2-1))-CX2(t2));
end

ystat1=abs((((1-t/N).*t/N).^delta).*(M1./t-M2./tflip));
[~,k]=max(ystat1);

if nargout==3
    Sigma=sqrt((MXX-M1(k)^2/k-M2(k)^2/(N-k))/(N-1))/sqrt(N);
end


% Sigma=std([(Xin(1:k)-mean(Xin(1:k))) (Xin(k+1:N)-mean(Xin(k+1:N)))])/sqrt(N);
% Sigma1=std([(X(t1:t1+k-1)-mean(X(t1:t1+k-1))) (X(t1+k:t1+N-1)-mean(X(t1+k:t1+N-1)))])/sqrt(N);
% Let consider var instead of std
% Sigma2=( mean([(X(t1:t1+k-1)-mean(X(t1:t1+k-1))) (X(t1+k:t1+N-1)-mean(X(t1+k:t1+N-1)))].^2) )/sqrt(N);
% c1=mean(X(t1:t1+k-1));
% c2=mean(X(t1+k:t1+N-1));
% N1=k;
% N2=N-k;
% Sigma3=mean( ([X(t1:t1+k-1) X(t1+k:t1+N-1)]-[c1*ones(1,N1) c2*ones(1,N2)]).^2)/sqrt(N);
% x=[X(t1:t1+k-1) X(t1+k:t1+N-1)];y=[c1*ones(1,N1) c2*ones(1,N2)];
% Sigma4= mean( (x-y).^2 )/sqrt(N);
% Sigma5 = sum (x.^2-2*x.*y +y.^2) /N/sqrt(N);
% Sigma6 = (sum(x.^2)/N - 2*sum(x.*y)/N + sum(y.^2)/N )/sqrt(N);
% Sigma7=( sum(X(t1:t1+N-1).^2) -N1*c1^2 -N2*c2^2 )/N/sqrt(N);
% c1=M1/N1,c2=M2/N2, MXX=sum(X(t1:t1+N-1).^2)
% Sigma8=((MXX-M1(k)^2/k-M2(k)^2/(N-k))/(N-1)) /sqrt(N);
% let divided by N instead of N-1, 
% Sigma9=sqrt(Sigma8)


function [tresh]=KSinverse(pf,M)
%[tresh]=KSinverse(pf,M)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   pf      rejection probability
%   M       number of terms in series
%
%   ystat   test statistic
%   k       estimated change point
%
%Estimate threshold for change points using inverse of Kolmogorov-Smirnov
%distribution

tresh=fzero(@KSdist,1,[],pf,M);

function [fun]=KSdist(Y,pf,M)
k=1:M;
s=(-ones(1,M)).^(k+1);
fun=2*sum(s.*exp(-(2*Y^2)*(k.^2)))-pf;


function Kout=finalal(CX1,CX2,CXX,Kin,m,delta)
%[meanK, Kout]=finalal(X,Xin,Kin,m,delta)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   X       Input data
%   Xin     Original time-series
%   Kin     vector containing estimated global change points
%   m       minimum interval size
%   delta   detect parameter
%
%   meanK   mean of values in intervals
%   Kout       final estimated change points
%
%Calculate final change point estimates and mean of each interval

if(nargin<3), m=5;  delta=0.5; end
if(nargin<4), delta=0.5; end
N=length(CX1); Z=length(Kin);

Kout=[];
if(length(Kin)>1)
    b=floor((Kin(2)+Kin(1))/2);
%     X1=X(1:b);
    [~,k1] = ystat(CX1,CX2,CXX,1,b,delta,m);

    Kout=[Kout; k1];
    for i=2:length(Kin)-1
        a=floor((Kin(i)+Kin(i-1))/2)+1;
        b=floor((Kin(i+1)+Kin(i))/2);
%         Xtemp=X(a:b);
        [~,ktemp] = ystat(CX1,CX2,CXX,a,b,delta,m);
        Kout=[Kout; a+ktemp-1];
    end
    a=round((Kin(Z)+Kin(Z-1))/2);
%     X1=X(a:N);
    [~,k1] = ystat(CX1,CX2,CXX,a,N,delta,m);
    Kout=[Kout; a+k1-1];
else [~,k1] = ystat(CX1,CX2,CXX,1,N,delta,m); 
    Kout=[Kout; k1];
end

function [Kout]=diagn(X,CX1,CX2,CXX,Kin,delta,m,thresh)
%[Knew,Kalt]=diagn(X,K,delta,m,thresh)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   X       Input data
%   K       vector containing estimated global change points
%   delta   detect parameter
%   m       minimum interval size
%   thresh  minimum statistic for change point characterization
%
%   Kout    output containing estimated global change points
%
%check estimated global change points for errors and perform update

%Initialize
N=length(CX1); Z=length(Kin); Kout=[];


if(length(Kin)>1)   %If inital change points found

    %Check first change point
    b=floor((Kin(2)+Kin(1))/2);
    [stat1,k1,Sigma] = ystat(CX1,CX2,CXX,1,b,delta,m); 
%     Sigma=std(X(1:b))/sqrt(b); 
    Level=thresh*Sigma*2;
    if(stat1(k1)>Level), Kout=[Kout; k1]; end

    %Check middle change points
    for i=2:length(Kin)-1
        a=floor((Kin(i)+Kin(i-1))/2)+1;
        b=floor((Kin(i+1)+Kin(i))/2);
        [stat,ktemp,Sigma] = ystat(CX1,CX2,CXX,a,b,delta,m); 
%         Sigma=std(X(a:b))/sqrt(b-a+1); 
        Level=thresh*Sigma*2;
        if(stat(ktemp)>Level), Kout=[Kout; a+ktemp-1]; end
    end

    %Check final change point
    a=round((Kin(Z)+Kin(Z-1))/2);
    [stat2,k2,Sigma] = ystat(CX1,CX2,CXX,a,N,delta,m); 
%     Sigma=std(X(a:N))/sqrt(N-a+1); 
    Level=thresh*Sigma*2;
    if(stat2(k2)>Level), Kout=[Kout; a+k2-1]; end

    
else                %If no change points yet found
    [stat1,k1,Sigma] = ystat(CX1,CX2,CXX,1,N,delta,m);
    %disp('test');
%     Sigma=std(X)/sqrt(N);
    if stat1(k1)>thresh*Sigma*2, Kout=[Kout; k1]; end
end

% PalarmF(CX1,CX2,CXX,t1,t2,Kin,P,delta,thresh,m,epsilon)
function Kout=diagn2(CX1,CX2,CXX,Kin,delta,m,thresh,epsilon,fac)
% Kout=diagn2(X,Kin,delta,m,thresh,epsilon,fac)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   X       Input data
%   K       vector containing estimated global change points
%   delta   detect parameter
%   m       minimum interval size
%   thresh  minimum statistic for change point characterization
%
%   Kout    output containing estimated global change points
%
% Check if there is any change poine in each plate
%Initialize
Z=length(Kin); Kout=Kin;

if ~isempty(Kin)   %If inital change points found
    %Check first change point
    b=Kin(1);%floor((Kin(2)+Kin(1))/2);
%     X1 = X(1:b);
    kout=PalarmF(CX1,CX2,CXX,1,b,[],1,delta,thresh,m,epsilon);
    Kout=[Kout; kout];
    %Check middle change points
    for i=2:length(Kin)-1
        a=Kin(i-1);%floor((Kin(i)+Kin(i-1))/2)+1;
        b=Kin(i);%floor((Kin(i+1)+Kin(i))/2);
%         Xi = X(a:b);
        kout=PalarmF(CX1,CX2,CXX,a,b,[],1,delta,thresh,m,epsilon);
        Kout=[Kout; a+kout-1];
    end

    %Check final change point
    a=Kin(Z);%round((Kin(Z)+Kin(Z-1))/2);
%     X2=X(a:end);
    kout=PalarmF(CX1,CX2,CXX,a,length(CX1),[],1,delta,thresh,m,epsilon);
    Kout=[Kout; a+kout-1];
end
Kout=sort(Kout);


function meanK=getmeanK(X,Kin)  
meanK=zeros(size(X));
meanK(1:Kin(1)-1)=mean(X(1:Kin(1)-1));

for i=2:length(Kin)        
    meanK(Kin(i-1):Kin(i)-1)=mean(X(Kin(i-1)+1:Kin(i)-1));
end
meanK(Kin(end):end) = mean(X(Kin(end):end));

    