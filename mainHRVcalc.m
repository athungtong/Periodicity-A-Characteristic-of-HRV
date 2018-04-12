function mainHRVcalc(filepath,fOutName,w)
%%
subjectpath = fileparts(filepath);
if ~isfilefound(filepath),return;end
load(filepath); %Will get Events from this line
if ~isfield(Events,'RR_interval'),return;end

if isfilefound(fullfile(subjectpath,fOutName))
    load(fullfile(subjectpath,fOutName));
else
    AllIndex=struct;
end

if isfield(AllIndex,'HRV')
    if isfield(AllIndex.HRV,'RR')
        AllIndex.HRV = rmfield(AllIndex.HRV,'RR');
    end
end

w=2.5;
flagR=0;
SDANN=[];
RMSSD=[];
TRI=[];

%Time domain analysis: SDANN, RMSSD, Triangular index
Seg=floor(Events.RR_time(end)/60/w); %number of section, round up
for s=0:Seg-1
    t_index=find(Events.RR_time>w*60*s & Events.RR_time<=w*60*(s+1));
    if isempty(t_index),flagR=s; continue; end
    RR_interval = Events.RR_interval(t_index)*1e3; %convert to millisecond
    SDANN(s+1,:)=std(RR_interval);
    dRR = diff(RR_interval);
    RMSSD(s+1,:)=sqrt(mean(dRR.^2)); %This measure high freq variation in HR, reflect para symphathetic nervous syterm 

    dx=1/256*1e3;
    minx = 0.1*1e3;
    maxx = 0.6*1e3;
    [f edges]=gethistogram(RR_interval,dx,minx,maxx);
    TRI(s+1,:)=sum(f)/max(f);

end
AllIndex.HRV.RR.SDANN=SDANN;
AllIndex.HRV.RR.RMSSD=RMSSD;
AllIndex.HRV.RR.TRI=TRI;

%Frequency domain analysis:LFP, HFP, LHR
w=2.5;
flagR=0;
LFP=[];
HFP=[];
LHR=[];
TF=[0.1 3.5];
LF=[0.1 1.0];
HF=[1.0 3.5];

Seg=floor(Events.RR_time(end)/60/w); %number of section, round up
for s=0:Seg-1
    t_index=find(Events.RR_time>w*60*s & Events.RR_time<=w*60*(s+1));
    if isempty(t_index),flagR=s; continue; end
    RR_interval = Events.RR_interval(t_index);
    RR_time = Events.RR_time(t_index);

    Flomb = FlombCalculation(RR_time,RR_interval,TF,HF,LF);
    LFP(s+1,:)=Flomb.LFP;
    HFP(s+1,:)=Flomb.HFP;
    LHR(s+1,:)=Flomb.LHR;
end

AllIndex.HRV.RR.LF=LFP;
AllIndex.HRV.RR.HF=HFP;
AllIndex.HRV.RR.LHR=LHR;


%Nonlinear measures: SD1, SD2, SDR, DFA.A1, DFA.A2, SampEn
SD1=[];
SD2=[];
SDR=[];
CE=[];
alpha1=[];
alpha2=[];
SampEn=[];
dfa=[];
w=2.5;
for s=0:Seg-1
    t_index=find(Events.RR_time>w*60*s & Events.RR_time<=w*60*(s+1));
    if isempty(t_index),flagR=s; continue; end
    RR_interval = Events.RR_interval(t_index)*1e3; %convert to millisecond

    Poincare=PoincarePlotCalculation(RR_interval,1);
    SD1(s+1,:)=Poincare.SD1; 
    SD2(s+1,:)=Poincare.SD2;
    SDR(s+1,:)=Poincare.SDRatio; 
    
    [~, tau] = time_delayed_MI(RR_interval,0:10);
    [m maxnorm] = fnn(RR_interval,'max',tau,5,0);

    [tem1]=Entropy(maxnorm,center(RR_interval),m,tau,.3,'max');
    if (isfinite(tem1)), SampEn(s+1,:) = tem1; end
    
     alpha1(s+1,:)= DFA(RR_interval,16);
     alpha2(s+1,:)= DFA(RR_interval,64);
   
end
AllIndex.HRV.RR.SD1=SD1;
AllIndex.HRV.RR.SD2=SD2;
AllIndex.HRV.RR.SDRatio=SDR;

AllIndex.HRV.RR.Alpha1=alpha1;
AllIndex.HRV.RR.Alpha2=alpha2;
AllIndex.HRV.RR.SampEn=SampEn;

save(fullfile(subjectpath,fOutName),'AllIndex');
display('done')
