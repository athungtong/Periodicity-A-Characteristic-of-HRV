%%Compute and save basic HRV

eventpath = fullfile('CVC Events');
load(fullfile(eventpath,'fileinfo4.mat'));
numGroup=[zeros(7,1);ones(8,1)]; %1 for BVV 0 for CMV
filename='Events';
fOutName='BasicIndex';
w=2.5;

for i=1:length(fileinfo.subjectlist)
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline',[filename '.mat']) 
    mainHRVcalc(filepath,[fOutName '.mat'],w);
    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',['Events' ventpart '.mat']);
        mainHRVcalc(filepath,[fOutName ventpart '.mat'],w);
    end
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','Events.mat');   
    mainHRVcalc(filepath,[fOutName '.mat'],w);
end
tempMergeIndex('HRV','RR');

%% work on stat
[HRVall bsl rec xtickvec xlabelvec indices]= loadindex('BasicIndex','HRV','RR',1,0);
HRVBVV=HRVall(HRVall(:,end-1)==0,:);
HRVCMV=HRVall(HRVall(:,end-1)==1,:);
dummytick = [{''},{''},{''},{''},{''},{''}];

%% work on ploting
figure(2)
[ha, pos] = tight_subplot(4,3,[.03 .03],[.03 .03],[.01 .01]);
% [left right],[buttom, above], [left right of all subplot]
yl = [0 20;0 20;0 10; 0 100; 0 100; 0 10; 0 2; 0 2; 0 2;0 1;0 2;0 5];
tk=[-5 -5 -2.5 -0.25 -0.25 -2.5 -5/1000 -5/1000 -0.5 -0.5 -0.5 -0.5]/2;
unit={'ms.','ms.','','%','%','','ms.','ms.','','','',''};
set(gcf, 'Renderer', 'Painters');
group=unique(HRVBVV(:,end));
fid=fopen('newHRVstat.txt','w');

for numindex=1:12   
    axes(ha(numindex))
    h=ha(numindex);
    %pos=get(h,'pos'); pos(4) = pos(4)+0.035;
%     switch numindex
%         case{1}
%             pos1=get(h,'pos');
%         case{2}
%             pos2=get(h,'pos');
%         case{3}
%             
%             pos(1)=pos1(1);
%         case{4,6}
%             pos(1)=pos2(1);
%         case{5}
%             pos(1)=pos1(1)+0.22;
%     end
    
%     set(h,'pos',pos)
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
    set(gca,'XTick',xtickvec+.11,'XTickLabel',dummytick,'fontsize',8)
    
%     set(gca,'XTick',xtickvec+.11,'XTickLabel',xlabelvec,'fontsize',7) 
    ylim(yl(numindex,:)); 
    title(indices{numindex}(1:end-3),'fontsize',12); 
    ylabel(unit(numindex),'fontsize',10);
    if(numindex<3)
        %title(['log(' indices{numindex}(1:end-3) ')'],'fontsize',10); 
    end


    xlim([group(1)-0.5 group(end)+0.5]);
%     if numindex~=1 && numindex~=6
%         set(gca,'yscale','log');
%         ylim([5e-4 5e-1]);
%         set(gca,'YTick',[1e-3 1e-2 1e-1]);
%     end
    for i=1:length(group)
        maxy = get(gca,'ylim');
%         text(group(i)+.11,maxy(2)/1.5,num2str(pvalue(i),3),'Rotation',45);
    end
%     hold on;
%     if numindex==1,ylim([0 0.801]);end;
%     if numindex==2,ylim([0 0.601]);end;
%     if numindex==4,ylim([0 0.6]);end;
%     ylim([0 1]); if numindex==5,ylim([0 3]);end
    yt=get(gca,'ytick'); maxy=yt(1);
    
    ypos = -0.15;
    for t=1:length(xlabelvec)
        %set(gca,'XTick',xtickvec+.11,'XTickLabel',xlabelvec)
         text(xtickvec(t),tk(numindex),xlabelvec{t},'rotation',45,'fontsize',8);
    end
    box off
    
%     if numindex==2 || numindex==3 ||numindex==4 ||numindex==5 %log scale
%         ytickvec=[log(1e-3) log(1e-2) log(1e-1)];
%         ylabelvec={'10^-3','10^-2','10^-1'};
%         set(gca,'YTick',ytickvec,'YTickLabel',ylabelvec)
%     end
end
fclose(fid);

fig1h=figure(2);
set(fig1h,'paperposition',[0.25 3 7.07 4.21]);
%print(fig1h,'basicHRV','-depsc');


