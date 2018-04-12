% function [SlopeEst,nTotal,FnTotal,ChangeInterval,FitData]= DFA(Y,TotalNumPoints)
function [SlopeEst,nTotal,FnTotal]= DFA(Y,TotalNumPoints)

% This function calulate the Detrended Fluctuation Analysis
% Ref: Peng C.K. "Quatifiction of scaling exponents and crossover phonomena
%      in nonstationary heartbeat time series"
%
% Inputs
% Y              : Input time series data
% TotalNumPoints : The default value is 40
%
% Outputs
% nToatl  : Calculated n values
% FnTotal : Calculated Fn Values

if size(Y,1)>size(Y,2)
    Y=Y';
end

if nargin<2
    TotalNumPoints=40;
end

Len = length(Y);

Y=cumsum(Y-mean(Y));
nTotal=fix(logspace(log10(10),log10(Len/10),TotalNumPoints));
% nTotal = 10:2:100;
nTotal(diff(nTotal)==0)=[];
FnTotal=[];

SumN=cumsum([1:max(nTotal)]); 
SumN2=cumsum( ( [1: max(nTotal) ]).^2);
ForwardSumY=cumsum(Y);
BackwardSumY=cumsum(fliplr(Y));
ForwardSumYN=cumsum(Y.*[1:Len]);
BackwardSumYN=cumsum(fliplr(Y).*[1:Len]);

nTotal(nTotal<2)=[];
    
for  n = nTotal;

    InvMatA=([SumN2(n) SumN(n);SumN(n) n]^(-1));

    if rem(Len,n)==0
        % there is no need for backward estimation
        Yn=zeros(1,Len);
        Para=InvMatA*[ForwardSumYN(n);ForwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = ForwardSumY(n*i)-ForwardSumY((i-1)*n) ;
            Para=InvMatA*[ForwardSumYN(i*n)-ForwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*(1:n)+Para(2);
        end
        FnTotal=[FnTotal sqrt(mean((Y-Yn).^2))];
    else
        % forward and backward estimation
        
        % Forward estimation
        Yn=zeros(1,fix(Len/n)*n);
        Para=InvMatA*[ForwardSumYN(n);ForwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = ForwardSumY(n*i)-ForwardSumY((i-1)*n) ;
            Para=InvMatA*[ForwardSumYN(i*n)-ForwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*[1:n]+Para(2);
        end
        
        TempFn=sqrt(mean((Y(1:fix(Len/n)*n)-Yn).^2))/2;
        
        
        % Backward estimation
        Yn=zeros(1,fix(Len/n)*n);
        Para=InvMatA*[BackwardSumYN(n);BackwardSumY(n)];
        Yn(1:n)=Para(1)*[1:n]+Para(2);

        for i=2:fix(Len/n)
            Temp = BackwardSumY(n*i)-BackwardSumY((i-1)*n) ;
            Para=InvMatA*[BackwardSumYN(i*n)-BackwardSumYN((i-1)*n)-(i-1)*n*Temp; Temp];
            Yn((1:n)+(i-1)*n)=Para(1)*[1:n]+Para(2);
        end
        TempFn2=sqrt(mean(( Y([rem(Len,n)+1:Len])-fliplr(Yn) ).^2))/2;
        
        FnTotal=[FnTotal TempFn2+TempFn];
        
    end

end



FnTotal=log10(FnTotal);
nTotal=log10(nTotal);
SlopeEst=detectslope(FnTotal,nTotal,0);


return;

Len=length(FnTotal);

Slopes =zeros(1,Len);

% Create a vector of pointwise slopes
for i=2:Len-1
    Slopes(i)=(FnTotal(i+1)-FnTotal(i-1))/(nTotal(i+1)-nTotal(i-1));
end
Slopes(1) = (FnTotal(2)-FnTotal(1))/(nTotal(2)-nTotal(1));
Slopes(Len) = (FnTotal(Len)-FnTotal(Len-1))/(nTotal(Len)-nTotal(Len-1));


subplot(211)
plot(nTotal,Slopes,'*')
hold on



Slopes(find(isnan(Slopes)))=[];

[Out,Kout]=Palarm(Slopes);

PlatStarts=[1 Kout];

PlatEnds=[Kout-1 Len];



FitData=[];
for i = 1:length(PlatStarts)

    X = nTotal(PlatStarts(i):PlatEnds(i));
    Y = FnTotal(PlatStarts(i):PlatEnds(i));
    % Regression
    a(1,1)=length(X);
    a(1,2)=sum(X);
    a(2,1)=a(1,2);
    a(2,2)=sum(X.*X);
    b(1)=sum(Y);
    b(2)=sum(X.*Y);
    c = a\b';
    Intercept(i)=c(1);
    SlopeEst(i)=c(2);
    FitData=[FitData SlopeEst(i)*X+Intercept(i)*ones(size(X))];

end

ChangeInterval=[PlatStarts;PlatEnds];

for i=1:length(PlatStarts)
    plot(nTotal(PlatStarts(i):PlatEnds(i)),Out(PlatStarts(i):PlatEnds(i)),'color','r','LineWidth',1.5)
    grid on
    hold on
end

xlabel('Log(n)','FontName','Arial','FontSize',12)
ylabel('Estimated Slopes','FontName','Arial','FontSize',12)
legend('Pointwise Slopes','Detected Slopes')
title('DFA Slopes','FontName','Arial','FontSize',14)




subplot(212)
plot(nTotal,FnTotal,'*')
hold on;
grid on;
plot(nTotal,FitData,'r')


function [Out,Kout]=Palarm(X,delta1,delta2,pf,m,M,epsilon)
%[Out,Kalt]=Palarm(X,delta1,delta2,pf,m,M,epsilon)
%
%   @2004 by Anatoly Zlotnik and Alexandra Piryatinska
%
%   Input   Default     Description
%   X                   Input data
%   delta1  0.8         detect parameter for phase I
%   delta2  0.5         detect parameter for phase II
%   pf      0.1        change point rejection probability
%   m       5           minimum interval size
%   M       1000        number of terms in KS series
%   epsilon 0.02        minimum relative distance of change point from 
%                       endpoints of local interval
%   
%   Out                 data with each interval replaced by local mean
%   Kout                vector containing estimated global change points
%  
%Estimation of change points and mean on each stable interval 
warning off MATLAB:fzero:UndeterminedSyntax

%scale input
Xin=X; Xbar=mean(X); X=(X-Xbar)/std(X);

if(nargin<2)
    delta1=0.8; delta2=0.5; pf=0.1; m=5; M=1000; epsilon=0.02;
end

if(size(X,1)>1) X=X'; end
if(size(X,1)>1) disp('error in input!'); return; end

%Calculate Initial Sigma estimate
N=length(X);
%X1=X(1:L*10);
%XX=reshape(X1,10,L);
%Sigma=mean(std(XX))/(L*10);

%Calculate Initial change point threshold
thresh=abs(KSinverse(pf,M));
%thresh=abs(KSinverse(pf,M))*Sigma;

%Calculate initial change points
P0=1; Kin0=[];
[Kone]=PalarmF(X,Kin0,P0,delta1,thresh,m,epsilon);
Kone=sort(Kone);
if(length(Kone)==0) Kout=Kone; Out=mean(Xin)*ones(N,1); return; end

%Calculate new Sigma estimate
%Y=X; Chng=[1 Kone N];
%for i=1:length(Chng)-1
%    Ytemp=Y(Chng(i):Chng(i+1)-1);
%    Y(Chng(i):Chng(i+1)-1)=Ytemp-mean(Ytemp);
%end
%Y(N)=Y(N)-Xbar;
%Sigma=std(Y)/sqrt(N);

%Calculate new threshold
thresh=abs(KSinverse(pf/10,M));

[Kout]=diagn(X,Kone,delta2,m,thresh);

[Out,Kout]=finalal(X,Xin,Kout,m,delta2,thresh);

if(size(Out,2)>1) Out=Out'; end


function [ystat1,k] = ystat(X,delta,m)
%[ystat1,k] = ystat(X,delta,m)
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

N=length(X); if(size(X,2)>1) X=X'; end
if(N<m) ystat1=zeros(1,N); k=1; return; end
M1=cumsum(X(1:N-1))./[1:N-1]'; 
M2=flipdim(cumsum(flipdim(X(2:N),1))./[1:N-1]',1);
ystat1=abs((((1-[1:N-1]'/N).*[1:N-1]'/N).^delta).*(M1-M2));
[M,k]=max(ystat1);



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
k=[1:M];
s=(-ones(1,M)).^(k+1);
fun=2*sum(s.*exp(-(2*Y^2)*(k.^2)))-pf;




function [Kout]=PalarmF(X,Kin,P,delta,thresh,m,epsilon)
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

[stat,k]=ystat(X,delta,m);
Sigma=std(X)/sqrt(length(X)); Level=thresh*Sigma;
D=max([round(epsilon*length(X)) 1]);
if(stat(k)<Level) Kout=Kin; return;
elseif(stat(k)>Level)
    X1=X(1:k-D);P1=P;
    X2=X(k+D:length(X));P2=P+k+D-1;
    if(length(X1)<=m) Ktemp1=Kin; elseif(length(X1)>m)
    Ktemp1=PalarmF(X1,Kin,P1,delta,thresh,m,epsilon); end
    if(length(X2)<=m) Ktemp2=Ktemp1; elseif(length(X2)>m)
    Ktemp2=PalarmF(X2,Ktemp1,P2,delta,thresh,m,epsilon); end
    if(length(X1)>m&length(X2)>m) Kout=[Ktemp2 P+k-1]; else
        Kout=Ktemp2; end
end



function [Kout]=diagn(X,Kin,delta,m,thresh)
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
N=length(X); Z=length(Kin); Kout=[];


if(length(Kin)>1)   %If inital change points found

    %Check first change point
    b=floor((Kin(2)+Kin(1))/2);
    [stat1,k1] = ystat(X(1:b),delta,m); Sigma=std(X(1:b))/sqrt(b); Level=thresh*Sigma;
    if(stat1(k1)>Level) Kout=[Kout k1]; end

    %Check middle change points
    for i=2:length(Kin)-1
        a=floor((Kin(i)+Kin(i-1))/2)+1;
        b=floor((Kin(i+1)+Kin(i))/2);
        [stat,ktemp] = ystat(X(a:b),delta,m); Sigma=std(X(a:b))/sqrt(b-a); Level=thresh*Sigma;
        if(stat(ktemp)>Level) Kout=[Kout a+ktemp-1]; end
    end

    %Check final change point
    a=round((Kin(Z)+Kin(Z-1))/2);
    [stat2,k2] = ystat(X(a:N),delta,m); Sigma=std(X(a:N))/sqrt(N-a); Level=thresh*Sigma;
    if(stat2(k2)>Level) Kout=[Kout a+k2-1]; end

    
else                %If no change points yet found
    [stat1,k1] = ystat(X,delta,m);
    %disp('test');
    if(stat1(k1)>thresh*std(X)/sqrt(N)) Kout=[Kout k1]; end
end
   

function [meanK, Kout]=finalal(X,Xin,Kin,m,delta,thresh)
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

if(size(X,2)>1), X=X'; end
if(nargin<3), m=5;  delta=0.5; end
if(nargin<4) delta=0.5; end
N=length(X); Z=length(Kin);

Kout=[];
if(length(Kin)>1)
    b=floor((Kin(2)+Kin(1))/2);
    X1=X(1:b);
    [stat1,k1] = ystat(X1,delta,m);
    if(k1>1)
        me=mean(Xin(1:k1-1));
        meanK=me*(ones(1,k1-1));
    else
        meanK=[];
    end
    Kout=[Kout k1];
    for i=2:length(Kin)-1
        a=floor((Kin(i)+Kin(i-1))/2)+1;
        b=floor((Kin(i+1)+Kin(i))/2);
        Xtemp=X(a:b);
        [stat,ktemp] = ystat(Xtemp,delta,m);
        Kout=[Kout a+ktemp-1];
        M1=mean(Xin(Kout(i-1)+1:Kout(i)-1));
        M2=M1*(ones(1,Kout(i)-Kout(i-1)));
        meanK=[meanK M2];
    end
    a=round((Kin(Z)+Kin(Z-1))/2);
    X1=X(a:N);
    [stat,k1] = ystat(X1,delta,m);
    Kout=[Kout a+k1-1];
    M1=mean(Xin(1+Kout(Z-1):Kout(Z)-1));
    M2=M1*(ones(1 ,Kout(Z)-Kout(Z-1)));
    meanK=[meanK M2];
    M3=mean(Xin(Kout(Z)+1:N));
    meanK=[meanK  M3*ones(1 ,N-Kout(Z)+1)];
else [stat1,k1] = ystat(X,delta,m); Kout=[Kout k1];
    meanK=[ones(1,k1)*mean(Xin(1:k1)) ones(1,N-k1)*mean(Xin(k1:N))]; 
end