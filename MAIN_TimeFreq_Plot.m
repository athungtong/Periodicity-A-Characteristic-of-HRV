                
clear all; 
eventpath = fullfile('CVC Events');

load(fullfile(eventpath,'fileinfo.mat'));
% w=200; %200 breaths
% w=2.5; % 2.5 minutes
tau=1;
fOutName='AllIndex'; couptype = 'HRV';  isplot=1;

% fOutName='Events'; %special for recording nm for each breath

technique{1}='TPVA'; w=100; nump=4; pMA=[1 10 20]; Maxtau=100; tau = 1;
freq = (0.1:1e-3:2); 
freq=fliplr(freq);
width = 3.475; 
%width = 7.134; % width of both column

height = 2.5;    % Height in inches
width=5;
alw = 0.75;    % AxesLineWidth
fsz = 8;      % Fontsize
lw = 0.72;      % LineWidth
msz = 6;       % MarkerSize

close all
figure(1)
% [ha, pos] = tight_subplot(3,2,[.03 .02],[.1 .05],[.08 .02]); % [left right],[buttom, above], [left right of all subplot]
[ha, pos] = tight_subplot(3,2,[.03 .06],[.12 .05],[.12 .02]); % [buttom, above],[buttom, above], [left right of all subplot]
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*100, height*100]); %<- Set size

tend=0;




flist=[8,12];
% for i=1:length(fileinfo.subjectlist)
% for k=1:length(flist)
    i=flist(1);
    HR=[];FR=[];VR=[];
    VT=zeros(1,fileinfo.numvent(i));
    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];
        display([fileinfo.Group(i,:) fileinfo.Subject(i,:) ventpart]);
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
        % tend=mainFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,tend,ha(1),ha(3),ha(5),[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart]);

        temp=size(VR,2);
        [HR FR VR]=mainTimeFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,VT,isplot,[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart],HR,FR,VR,freq);
        VT(j)=size(VR,2)-temp;
    end 

    VT=cumsum(VT);
    
    yel=colormap(flipud(gray));
    yel(:,1)=linspace(1,173/255,length(yel));
    yel(:,2)=linspace(1,118/255,length(yel));
    red=colormap(flipud(gray));
    red(:,1)=linspace(1,255/255,length(red));
    
    blu=colormap(flipud(gray));
    blu(:,2)=linspace(1,99/255,length(blu));
    blu(:,3)=linspace(1,255/255,length(blu));
    colormap(yel);
    
    axes(ha(1));    imagesc(1:size(VR,2)*2.5,freq,VR);
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    ylabel({'Ventilator','(Hz)'},'fontsize',8); 
    set(gca,'YDir','normal')
    set(gca,'XTickLabel','');
    set(gca,'YTick',0:0.5:2,'YTickLabel',{'0','','1','','2'},'fontsize',5);ylim([0,2]);
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    box off;

    axes(ha(3));imagesc(1:size(HR,2)*2.5,freq,HR);set(gca,'YDir','normal')
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    ylabel({'RR Interval','(Hz)'},'fontsize',8); 

    set(gca,'XTickLabel','');
    set(gca,'YTick',0:0.5:2,'YTickLabel',{'0','','1','','2'},'fontsize',5);ylim([0,2]);
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    box off;
    axes(ha(5));imagesc(1:size(FR,2)*2.5,freq,FR); set(gca,'YDir','normal'); xlabel('Time (minute)','fontsize',8);
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    set(gca,'YTick',0:0.5:2,'YTickLabel',{'0','','1','','2'},'fontsize',5);ylim([0,2]);
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    ylabel({'Breathing','(Hz)'},'fontsize',8); 
box off;
%     con=input('con'); if ~isempty(con),break;end
    
    
 %%%%%%%%%%%%%%%%%%%%%CMV

    i=flist(2);
    HR=[];FR=[];VR=[];
    VT=zeros(1,fileinfo.numvent(i));
    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];
        display([fileinfo.Group(i,:) fileinfo.Subject(i,:) ventpart]);
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
        % tend=mainFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,tend,ha(1),ha(3),ha(5),[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart]);

        temp=size(VR,2);
        [HR FR VR]=mainTimeFreqmod(filepath,[fOutName ventpart '.mat'],technique,w,tau,couptype,VT,isplot,[fileinfo.Group(i,:) ' ' fileinfo.Subject(i,:) ' ' ventpart],HR,FR,VR,freq);
        VT(j)=size(VR,2)-temp;
    end 

    VT=cumsum(VT);
    
    axes(ha(2));    imagesc(1:size(VR,2)*2.5,freq,VR);
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    %ylabel({'HR','(Beats/min)'},'fontsize',8); 
    set(gca,'YDir','normal')
    set(gca,'XTickLabel','');
    set(gca,'YTick',0:0.5:2,'YTickLabel','','fontsize',5)
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    box off;
    
    axes(ha(4));imagesc(1:size(HR,2)*2.5,freq,HR);set(gca,'YDir','normal')
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    %ylabel({'FR','(Brths/min)'},'fontsize',8); 
    set(gca,'XTickLabel','');
    set(gca,'YTick',0:0.5:2,'YTickLabel','','fontsize',5)
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    box off;   
    
    axes(ha(6));imagesc(1:size(FR,2)*2.5,freq,FR); set(gca,'YDir','normal'); xlabel('Time (minute)','fontsize',8);
    %hold on;for(j=1:fileinfo.numvent(i)),  plot([VT(j) VT(j)]*2.5,[0 2],'r'), end, hold off;
    %ylabel({'VR','(Cycles/min)'},'fontsize',8); 
%     set(gca,'YTick',0:0.5:2,'YTickLabel','','fontsize',5)
    set(gca,'YTick',0:0.5:2,'YTickLabel','','fontsize',5)
    set(gca,'XTick',VT*2.5); %grid on; set(gca,'YGrid','off')
    box off;
% end
fig1h=figure(1);
set(fig1h,'paperposition',[0.25 2.5 width height]);
print(fig1h,'TimeFreqyel','-depsc');


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