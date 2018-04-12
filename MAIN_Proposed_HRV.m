%%Compute and save basic HRV
eventpath = fullfile('CVC Events');
load(fullfile(eventpath,'fileinfo4.mat'));
numGroup=[zeros(7,1);ones(8,1)]; %1 for BVV 0 for CMV
filename='Events';
fOutName='AllIndex';
w=2.5;
isplot=0;
rawfile=[];
pairname='RR';
nump=4;
pMA = [1 10 100 500];
Maxtau = 100; tau=100;
label=[];
figpath='';
for i=1:length(fileinfo.subjectlist)
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline',[filename '.mat']) 
    mainAutocorr(filepath,rawfile,[fOutName '.mat'],pairname,nump,pMA,Maxtau,tau,w,isplot,label,figpath);

    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
        mainAutocorr(filepath,rawfile,[fOutName ventpart '.mat'],pairname,nump,pMA,Maxtau,tau,w,isplot,label,figpath);
    end
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','Events.mat');   
    mainAutocorr(filepath,rawfile,[fOutName '.mat'],pairname,nump,pMA,Maxtau,tau,w,isplot,label,figpath);
end
tempMergeIndex('HRV','RR');

%% work on stat
[HRVall bsl rec xtickvec xlabelvec indices]= loadindex('AllIndex','HRV','RR',1,0);
HRVBVV=HRVall(HRVall(:,end-1)==0,:);
HRVCMV=HRVall(HRVall(:,end-1)==1,:);
dummytick = [{''},{''},{''},{''},{''},{''}];

%% work on ploting
figure(2)
[ha, pos] = tight_subplot(2,3,[.03 .03],[.03 .03],[.01 .01]);
% [left right],[buttom, above], [left right of all subplot]
yl = [0 20;0 20;0 10; 0 100; 0 100; 0 10; 0 2; 0 2; 0 2;0 1;0 2;0 5];
tk=[-5 -5 -2.5 -0.25 -0.25 -2.5 -5/1000 -5/1000 -0.5 -0.5 -0.5 -0.5]/2;
%unit={'ms.','ms.','','%','%','','ms.','ms.','','','',''};
set(gcf, 'Renderer', 'Painters');
group=unique(HRVBVV(:,end));
fid=fopen('newHRVstat.txt','w');

for numindex=1:5   
    axes(ha(numindex))
    h=ha(numindex);
    x=[HRVBVV(:,numindex) HRVBVV(:,end) HRVBVV(:,end-4)]; 
    x(x(:,1)==Inf,:)=[]; %x(:,1)=sqrt(x(:,1));
    x(isnan(x),:)=[];
    xi = cell(length(group),1);
    for i = 1:length(group)
        xi{i} = x(x(:,2)== group(i),[1 3]);
    end
    
    y=[HRVCMV(:,numindex) HRVCMV(:,end)  HRVCMV(:,end-4)]; 
    y(y(:,1)==Inf,:)=[]; %y(:,1)=sqrt(y(:,1));
    y(isnan(y),:)=[];
    yi = cell(length(group),1);
    for i = 1:length(group)
        yi{i} = y(y(:,2)== group(i),[1 3]);
    end
    pvalue = zeros(length(group),1);
    fprintf(fid,'%s\r\n',['Index: ' indices{numindex}(1:end-3)]);
    %Try tramsformation
    if(numindex==4 || numindex==5)
        for i=1:length(group) 
            xi{i}=100*(xi{i});
            yi{i}=100*(yi{i});
        end
    end
    for i = 1:length(group)
        pvalue(i) = mainmanova(xi{i},yi{i},h,group(i),fid,char(xlabelvec(i)));
    end
    %set(gca,'XTick',xtickvec+.11,'XTickLabel',dummytick,'fontsize',8)
    %ylim(yl(numindex,:)); 
    title(indices{numindex}(1:end-3),'fontsize',12); 
    ylabel(unit(numindex),'fontsize',10);

    xlim([group(1)-0.5 group(end)+0.5]);
    for i=1:length(group)
        maxy = get(gca,'ylim');
    end
    yt=get(gca,'ytick'); maxy=yt(1);
    
    ypos = -0.15;
    for t=1:length(xlabelvec)
         text(xtickvec(t),tk(numindex),xlabelvec{t},'rotation',45,'fontsize',8);
    end
    box off
    
end
fclose(fid);

fig1h=figure(2);
set(fig1h,'paperposition',[0.25 3 7.07 4.21]);
%print(fig1h,'basicHRV','-depsc');


