function [Flomb f A]=FlombCalculation(RR_time,RR_interval,TF,HF,LF)
%%
% input 
% RR_time =time when beat2beat occur
% RR_interval=beat to beat interval
% TF: total freq range
% HF high freq range
% LF low freq range
if nargin<3
    TF = [0.04 0.4];
end
if nargin<4
    LF=[min(TF) mean(TF)];
    HF=[mean(TF) max(TF)];
% HF=[.15 .4]; %High Freq range
% LF=[.04 .15];%Low Freq range. For Human
% HF=[.25 1.5]; %High Freq range
% LF=[.01 .1];%Low Freq range. 
end
X=(RR_time);
Y=(RR_interval);
X=X-min(X);
[WK1 WK2]=FASPER(X,Y);
f=WK1;A=WK2;

% plot(f,(A/sum(A)));
% con=input('con?');

LFsection = A((f>LF(1)) & ((f<LF(2))));
HFsection = A((f>HF(1)) & ((f<HF(2))));
TFsection = A((f>=TF(1)) & ((f<TF(2))));
if ~isempty(LFsection)
%     LFP=max(LFsection)/(sum(A)); 
    LFP=sum(LFsection)/sum(A);
%     LFP=sum(LFsection)/sum(A)/range(A);
else
    LFP = -Inf;
end

if ~isempty(HFsection)
%     HFP=max(HFsection)/(sum(A)); 
    HFP=sum(HFsection)/sum(A);
%     HFP=sum(HFsection)/sum(A)/range(A);
else
    HFP = -Inf;
end

if ~isempty(TFsection)
%     TFP=max(TFsection)/(sum(A)); 
    TFP=sum(TFsection)/sum(A);
%     TFP=sum(TFsection)/sum(A)/range(A);
else
    TFP = -Inf;
end


LHR=LFP./HFP;
Flomb.LFP=LFP;
Flomb.HFP=HFP;
Flomb.TFP=TFP;
Flomb.LHR=LHR;

end

