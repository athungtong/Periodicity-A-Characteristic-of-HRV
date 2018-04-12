                
% clear all; 
close all
eventpath = fullfile('CVC Events');
figure(1)
height = 2.5;    % Height in inches
width=5;
alw = 0.75;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 0.72;      % LineWidth
msz = 6;       % MarkerSize

% [ha, pos] = tight_subplot(3,2,[.03 .02],[.1 .05],[.08 .02]); % [left right],[buttom, above], [left right of all subplot]
[ha, pos] = tight_subplot(3,2,[.03 .06],[.12 .05],[.12 .02]); % [buttom, above],[buttom, above], [left right of all subplot]
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size


load(fullfile(eventpath,'fileinfo.mat'));
% w=200; %200 breaths
% w=2.5; % 2.5 minutes
tau=1;
fOutName='AllIndex'; couptype = 'HRV';  isplot=1;

% fOutName='Events'; %special for recording nm for each breath

technique{1}='TPVA'; w=100; nump=4; pMA=[1 10 20]; Maxtau=100; tau = 1;
freq = (0.1:1e-3:2); 
freq=fliplr(freq);
%fig=figure('Position', [100, 100, 1200, 700]);
tend=0;
colormap(flipud(gray));

i=8;
HR=[];FR=[];VR=[];
VT=zeros(1,fileinfo.numvent(i));
j=4;
ventpart=['v' char(num2str(j))];
display([fileinfo.Group(i,:) fileinfo.Subject(i,:) ventpart]);
filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
tend=mainFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,tend,ha(1),ha(3),ha(5),[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart]);

axes(ha(1))
set(gca,'XTick',0:0.5:2,'XTickLabel','','fontsize',5);xlim([0,2]);
set(gca,'YTick',0:0.05:0.1,'YTickLabel',{'0','0.05','0.1'},'fontsize',5);ylim([0,0.1]);
ylabel('Incidence','fontsize',8);

axes(ha(3))
set(gca,'XTick',0:0.5:2,'XTickLabel',{'0','','1','','2'},'fontsize',5);xlim([0,2]);
set(gca,'YTick',0:0.025:0.05,'YTickLabel',{'0','0.025','0.05'},'fontsize',5);ylim([0,0.05]);
ylabel('nrmlzd PSD','fontsize',8); 
xlabel('Frequency (Hz)','fontsize',8);

axes(ha(5))
set(gca,'XTick',0:25:100,'XTickLabel',{'0','25','50','75','100'},'fontsize',5);xlim([0,100]);
set(gca,'YTick',-1:1:1,'YTickLabel',{'-1','0','1'},'fontsize',5);ylim([-1 1]);
ylabel('r','fontsize',8); xlabel('Delay (samples)','fontsize',8);


i=12;
HR=[];FR=[];VR=[];
VT=zeros(1,fileinfo.numvent(i));
j=5;
ventpart=['v' char(num2str(j))];
display([fileinfo.Group(i,:) fileinfo.Subject(i,:) ventpart]);
filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
tend=mainFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,tend,ha(2),ha(4),ha(6),[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart]);

axes(ha(2))
set(gca,'XTick',0:0.5:2,'XTickLabel','','fontsize',5);xlim([0,2]);
set(gca,'YTick',0:0.5:1,'YTickLabel',{'0','0.5','1'},'fontsize',5);ylim([0,1]);

axes(ha(4))
set(gca,'XTick',0:0.5:2,'XTickLabel',{'0','','1','','2'},'fontsize',5);xlim([0,2]);
set(gca,'YTick',0:0.25:0.5,'YTickLabel',{'0','0.25','0.5'},'fontsize',5);ylim([0,0.5]);
xlabel('Frequency (Hz)','fontsize',8);

axes(ha(6))
set(gca,'XTick',0:25:100,'XTickLabel',{'0','25','50','75','100'},'fontsize',5);xlim([0,100]);
set(gca,'YTick',-1:1:1,'YTickLabel',{'-1','0','1'},'fontsize',5);ylim([-1 1]);
xlabel('Delay (samples)','fontsize',8);

fig1h=figure(1);
set(fig1h,'paperposition',[0.25 2.5 width height]);
%print(fig1h,'Histogram','-depsc');

%%
fOutName='AllIndex1';
% fOutName='Events1';
parfor i=17:length(fileinfo.subjectlist) 
    tend=0;
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','Events1.mat');    
    rawfile=fullfile(datapath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','baseline1.mat');    
%     tend=RVCcalc(filepath,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
%     tend=mainPhaseSync(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
   tend=mainCoupDirection(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);

    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','Events1.mat');    
    rawfile=fullfile(datapath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','recovery1.mat');
%     tend=RVCcalc(filepath,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
    tend=mainCoupDirection(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
%     tend=mainPhaseSync(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);

end

fOutName='AllIndex2';
% fOutName='Events2';
parfor i=17:length(fileinfo.subjectlist)  
    tend=0;
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','Events2.mat');    
    rawfile=fullfile(datapath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','baseline2.mat');
%    tend=RVCcalc(filepath,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
%     tend=mainPhaseSync(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
      tend=mainCoupDirection(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);

    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','Events2.mat');    
    rawfile=fullfile(datapath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','recovery2.mat');
%    tend=RVCcalc(filepath,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
   tend=mainCoupDirection(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);
%       tend=mainPhaseSync(filepath,rawfile,[fOutName '.mat'],technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau);

end
%%
pairname=cell(4,1);
if strcmp(couptype,'CRC')
    pairname{1}='RR';pairname{2}='RI';pairname{3}='ER';pairname{4}='RE'; numindex=3;
elseif strcmp(couptype,'RVC')
    pairname{1}='UI';pairname{2}='IU';pairname{3}='DE';pairname{4}='ED'; numindex=4;
elseif strcmp(couptype,'CVC')
    pairname{1}='UR';pairname{2}='DR'; numindex = 4;
end

for i=1:length(pairname)  
  if any(strcmp(technique,'TPVA'))
      tempMergeTPVAIndex(pairname{i},'TPVA',length(pMA),nump,1)
      tempMergeTPVAIndex(pairname{i},'TPVTD',length(pMA),nump,Maxtau)
  else
      tempMergeIndex(couptype,pairname{i});
  end 
end

%%
fOutName='AllIndex';couptype='RVC';
ismultivent=1; normalize=0;

for i=1:length(pairname) 
    if any(strcmp(technique,'TPVA'))
        mainTPVAstat(pairname{i},(pMA),nump,1,ismultivent,i);
        mainTPVTDstat(pairname{i},4,Maxtau,ismultivent,i);
    else
        mainCRCstat(fOutName,couptype,pairname{i},ismultivent,normalize,numindex,i);
    end

end