function tend=mainCoupDirection(filepath,rawfile,fOutName,technique,w,tau,couptype,tend,isplot,nump,pMA,Maxtau)

tend=0;

if ~isfilefound(filepath),return;end
subjectpath = fileparts(filepath);
load(filepath);

x='U';y='I'; 
% x='D';y='E'; nx=1:7; my=1:7;

if ~isfield(Events,[x x '_time']),return;end 
if ~isfield(Events,[y y '_time']),return;end

PlatStarts=Events.UInmInfo.PlatStarts;
PlatEnds=Events.UInmInfo.PlatEnds;
optnm=Events.UInmInfo.optnm;
tvdat=Events.UInmInfo.tvdat;


tx=getfield(Events,[x x '_time']);
ty=getfield(Events,[y y '_time']);
w=2.0;
nbin = 16; tau=2:2:10; fignum=1;
safr=100; 
Ns=19;
phx=(0:length(tx)-1)*2*pi;
phy=(0:length(ty)-1)*2*pi;
[txy, phx, phy]=samplingphase(tx,ty,phx,phy,1/safr);    
Seg=floor(txy(end)/w/60);
if Seg==0,Seg=1;end
tcoup=60*w*(0.5:(Seg-0.5))';

D=-10*ones(Seg,1); 
noD=D;  uniD12=D;   uniD21=D;   biD=D; MI2=D;  nmr=D; hMI=D;
I12=D;  I21=D;  I12s=zeros(Ns,size(D,1)); I21s=I12s; MI = D; MIs = I12s; 

n=0; %use to count number of good segment
for i=1:length(PlatStarts)
    dt=tvdat(PlatEnds(i))-tvdat(PlatStarts(i)); %this is the total time in this plate
    if w>dt/60,continue;end %if w > window, then skip
    t0=tvdat(PlatStarts(i));
    Seg=floor(dt/w/60);
    for s=0:Seg-1      
        indextxy= txy-t0>=w*s*60 & txy-t0<w*(s+1)*60;

        indextx= tx-t0>=w*s*60 & tx-t0<w*(s+1)*60;
        indexty= ty-t0>=w*s*60 & ty-t0<w*(s+1)*60;
        if isempty(phx(indextxy)),continue; end

        [~, I12(n+1) I21(n+1) theta1 theta2 Dtheta1 Dtheta2 MI(n+1)]=main_coupdirecinfo(phx(indextxy),phy(indextxy),safr,tau,nbin);
        [~, I12s(:,n+1) I21s(:,n+1),MIs(:,n+1)]=direcsurrogate_time(tx(indextx),ty(indexty),safr,tau,nbin,Ns,fignum,theta1,theta2,Dtheta1,Dtheta2,0);
        [noD(n+1) uniD12(n+1) uniD21(n+1) biD(n+1) D(n+1) MI2(n+1) hMI(n+1)]=direcsurrogatetest(I12(n+1),I21(n+1),I12s(:,n+1),I21s(:,n+1),MI(n+1),MIs(:,n+1));    
        nmr(n+1)=optnm(i,1)/optnm(i,2);
        n=n+1;tcoup(n)=t0+60*w*(s+0.5); 
        
    end
end

toberem = n+1:length(D);

% toberem=find(D==-10);
tcoup(toberem)=[];
I12(toberem)=[];
I21(toberem)=[];
MI(toberem)=[];
MI2(toberem)=[];
hMI(toberem)=[];
uniD12(toberem)=[];
uniD21(toberem)=[];
biD(toberem)=[];
D(toberem)=[];
noD(toberem)=[];
nmr(toberem)=[];

I12s(:,toberem)=[];
I21s(:,toberem)=[];
MIs(:,toberem)=[];

summaryDirec=[length(noD); sum(noD); sum(uniD12); sum(uniD21); sum(biD)];
biDirecScore = D(biD==1);

dir = 1*noD+2*uniD12+3*uniD21+4*biD; %dir=1 for noD, 2 for uniD12

biscore=zeros(size(D));
for i=1:length(D)
    if biD(i)==1,biscore(i)=D(biD(i));else biscore(i)=0;end
end
summaryDirec2=[dir biscore hMI nmr];

if isplot
    xlabelvec={};
    for i=1:length(tcoup)
        xlabelvec{i}=num2str(tcoup(i)/60+w/2);
    end
    
    figure(11)
    subplot(3,2,1:2);
    plot(tx(1:end-1)/60,diff(tx),'c.-','markersize',6); hold on;
    plot(ty(1:end-1)/60,diff(ty),'g.-','markersize',6); 
    
    % add frequency ratio information
    for i=1:length(PlatStarts)
        S=[num2str(optnm(i,1)) ':' num2str(optnm(i,2))];
        for j=1:1
            subplot(3,1,j);
            text(tvdat(PlatStarts(i))/60,0.7,S,'fontsize',9,'rotation',90);    grid on; 
        end
    end  
    % add sync index based on cir var of vdat
    if isfilefound(fullfile(subjectpath,fOutName))
        load(fullfile(subjectpath,fOutName)); %This will load AllIndex   
        syncind=AllIndex.RVC2.IR.cvar1;
        tindex=60*1*(0.5:(length(syncind)-0.5));
        plot(tindex/60,syncind,'.-');
    end
    hold off; grid on;
    title('xx interval and yy interval'); legend('x','y');

    subplot(323)
      plot(tcoup/60,I12,'o-'); hold on; 
   boxplot(I12s,'position',tcoup/60);hold off; 
    grid on;  ylim auto;%ylim([0.0 2]);
    title('I12 and I12s'); 
    set(gca,'XTick',tcoup/60+w/2,'XTickLabel',xlabelvec)
    
    subplot(324)
    plot(tcoup/60,I21,'o-'); hold on; 
    boxplot(I21s,'position',tcoup/60); hold off; 
    grid on;  ylim auto;%ylim([0.0 2]);
    title('I21 and I21s'); 
    
 set(gca,'XTick',tcoup/60+w/2,'XTickLabel',xlabelvec)
    subplot(325)
    plot(tcoup/60,MI,'o-'); hold on; 
    boxplot(MIs,'position',tcoup/60); hold off; 
    grid on;  ylim auto;%ylim([0.0 2]);
    title('MI and MIs'); 
 set(gca,'XTick',tcoup/60+w/2,'XTickLabel',xlabelvec)
 
    dir = [noD'; uniD12'; uniD21'; biD'];
    mark={'o','>','<','diamond'}; c={'k','b','r','g'};
    subplot(326)
    for i = 1:length(tcoup)
        if biD(i)==1,
            plot(tcoup(i)/60,D(i),'g','marker','diamond'); hold on;
        else
            plot(tcoup(i)/60,0,'marker',mark{find(dir(:,i))},'color',c{find(dir(:,i))},'markersize',8); hold on;
        end
    end
    hold off;grid on;  ylim auto;%ylim([0.0 2]);
    title('D'); 
    
    input('continue?');
end


if isfilefound(fullfile(subjectpath,fOutName))
   load(fullfile(subjectpath,fOutName)); %This will load AllIndex
else
    AllIndex=struct;
end

if isfield(AllIndex,'DIR')
    if isfield(AllIndex.DIR,([x y]))
        AllIndex.DIR = rmfield(AllIndex.DIR,([x y]));
    end
end

AllIndex.DIR.([x y]).summaryDirec=summaryDirec;
AllIndex.DIR.([x y]).summaryDirec2=summaryDirec2;
AllIndex.DIR.([x y]).biDirecScore=biDirecScore;
AllIndex.DIR.([x y]).MI=[tcoup MI2];

if ~isplot
    save(fullfile(subjectpath,fOutName),'AllIndex');
end



