function mainAutocorr(eventfile,rawfile,fOutName,pairname,nump,pMA,Maxtau,tau,w,isplot,label,figpath)
%%
% c1=[153 51 0]/255;
% c2=[174 119 0]/255;
% c3=[222 125 0]/255;
% c4=[255 200 100]/255;

subjectpath = fileparts(eventfile);
if ~isfilefound(eventfile),return;end
load(eventfile); %Will get Events from this line
if ~isfield(Events,'RR_interval'),return;end


tr=Events.RR_time; 
rr=Events.RR_interval;

Seg=floor(Events.RR_time(end)/60/w); %number of section, round up
if Seg==0,Seg=1;end

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

flagR = -1;
meanACF=[];
maxACF=[];
SampEn=[];
maxPER=[];
tRSE=[];
for s=0:Seg-1
    indexrr=find(tr>=w*60*s & tr<w*60*(s+1));
    if isempty(indexrr),flagR=s; continue; end
    RR = rr(indexrr);

    %%%% Working with autocorrelation of RR %%%%%%%%
    
    newRR = mydetrend(RR,1,tr(indexrr));
    [gammar d]=acf(newRR,round(tau));
    gammar(1:tau+1)=[]; 
    meanACF(s+1,1) = mean(abs(gammar(25:50)));
    maxACF(s+1,1) = max(abs(gammar(25:50)));
    
    [~, t] = time_delayed_MI(gammar,0:10);
    [m maxnorm] = fnn(gammar,'max',t,5,0);
    SampEn(s+1,1)= Entropy(maxnorm,center(gammar),m,t,.4,'max');      

%%%%%%%%% Working with Periodogram of RR interval %%%%%%%%%%%%%%%%
        [f P]=FASPER(tr(indexrr),RR,1);
        f(P<0)=[];  P(P<0)=[]; 
        P(f<0.1)=0;
        relatP = P/sum(P); 
        maxPER(s+1,1) = max(relatP);

        bin=20;
        dx = range(f)/bin;
        edges = min(f):dx:max(f);   
        edges(end)=Inf;
        P2=zeros(bin,1);
        for i = 1:bin
            P2(i) = sum(relatP(f>=edges(i) & f<edges(i+1)));
        end
        edges(P2==0)=[];
        P2(P2==0)=[];
        tRSE(s+1,1) = 1+dot(P2,log(P2))/log(bin);
        

    if isplot
        Fi=1./Events.II_interval;
        figname=[label ' minute ' num2str((s+0)*w) '-' num2str((s+1)*w)];

%%%% Plot Example %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         flag=2;
%         figure(2);
%         subplot(3,4,0+flag),plot(tr(indexrr)/60,60./RR,'-','color',c4); axis 'tight'; hold on;
%         set(gca,'fontsize',5); 
%         plot(tr(indexrr)/60,60./RR,'.','color',c1,'markersize',4); axis 'tight'; hold off;
%         ylabel('InstHR','fontsize',8);
%         xlabel('time (min)','fontsize',8);  
%         title(figname,'fontsize',8);
%         
%         subplot(3,4,4+flag),plot(f,(relatP),'color',c3); 
%         set(gca,'fontsize',5);
%         ylabel({'Periodogram','of RR interval'},'fontsize',8);
%         xlabel('Frequency (Hz)','fontsize',8);
%         xlim([0 max(f)]);   ylim([0 max(relatP)]);
%         text(0.1,max(relatP)/1.25,['tRSE = ' num2str(tRSE_P(s+1,1),2)],'fontsize',6); 
%         
%         
%         subplot(3,4,8+flag),plot(d(tau+2:end),(gammar),'-','color',c4);hold on;
%         plot(d(tau+2:end),(gammar),'.','color',c1,'markersize',4);hold off;
%         set(gca,'fontsize',5); 
%         xlabel('delay (samples)','fontsize',8);
%         ylabel({'Autocorrelation','of RR interval'},'fontsize',8);
%         ylim([-1 1]); grid off; xlim([0 tau]);
%         text(tau/20,0.8,['max ACF = ' num2str(maxACF3(s+1,1),2)],'fontsize',6); 
%         con = input('continue?');
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(1);
%         subplot(4,2,[3]),plot(tr(indexrr)/60,60./RR,'-','color','c'); axis 'tight'; hold on;
%         set(gca,'fontsize',5);
%         plot(tr(indexrr)/60,60./RR,'.','color','b'); axis 'tight'; hold off;
%         ylabel('InstHR','fontsize',8);
%         xlabel('time (min)','fontsize',8);

        plotFlomb(Events,RR,indexrr,tr,label,s,w,f,P,Fi)
%         subplot(325)
%         subplot(325);  hold on; plot(f,P2,'r');hold off;
% %         title([num2str(max(P)) ', ' num2str(temp.DFA1) ', ' num2str(temp.DFA2) ', ' num2str(temp.DFA3) ', '],'fontsize',8);
        h1=subplot(3,2,[1:2]);
        
        plot(d(tau+2:end),(gammar),'.-','color',c2);
        set(h1,'fontsize',5); 
        xlabel('delay (samples)','fontsize',8);
        ylabel({'Autocorrelation','of RR interval'},'fontsize',8);
        ylim([-1 1]); grid on;
        mn = mean(abs(gammar)); st=std(abs(gammar)); cv = st/mn;
%         title([' abs mn=' num2str(mn) ',abs st=' num2str(st) ',cv=' num2str(cv)],'fontsize',8);
        title(figname,'fontsize',8);
        subplot(3,2,[4 6]);
        
        plot3(gammar(1:end-2*t),gammar(1+t:end-t),gammar(1+2*t:end),'-','color',c4,'linewidth',0.05);hold on;
        plot3(gammar(1:end-2*t),gammar(1+t:end-t),gammar(1+2*t:end),'.','markersize',4,'color',c1,'linewidth',0.05);hold off;
        set(gca,'box','on','CameraPosition', [8.6097   10.2727    8.1791],'fontsize',5);
        zlabel({'TimeDelay Reconst','3D w/ unit delay'},'fontsize',8);
        set(gca,'xlim',[-1 1],'ylim',[-1 1],'zlim',[-1 1]);
        

%         fig1h=figure(1);
%         set(fig1h,'paperposition',[0.25 2.5 6 4]);
%         print(fig1h,'allplotCMV','-depsc');

        % tmir = (0:round(tau/sp))*sp;
        % subplot(326),plot( tmir(2:end),MIr(2:end),'.-');
        con = input('continue?');
    end
end


tRSE(1:flagR+1)=[];AllIndex.HRV.RR.tRSE=tRSE;
maxPER(1:flagR+1)=[];AllIndex.HRV.RR.maxPER=maxPER;

maxACF(1:flagR+1)=[];AllIndex.HRV.RR.maxACF=maxACF;
meanACF(1:flagR+1)=[];AllIndex.HRV.RR.meanACF=meanACF;

SampEn(1:flagR+1)=[];AllIndex.HRV.RR.SampEn=SampEn;

if ~isplot
    save(fullfile(subjectpath,fOutName),'AllIndex');
end

end

function plotFlomb(Events,RR,indexrr,tr,label,s,w,f,P,Fi)
    
    c1=[153 51 0]/255;  c2=[174 119 0]/255; c3=[222 125 0]/255;
    clear h1 h2;
    figname=[label ' minute ' num2str((s+0)*w) '-' num2str((s+1)*w)];
    figure(1);
    h{2}=subplot(325);
    set(h{2},'fontsize',5);
    totalp = sum(P);
    relatP = P/sum(P); 
    relatP(f<0.1)=0; 
%     P(f<0.1)=0;
    plot(f,(relatP),'color',c2); 
    title({'Periodogram of RR interval'},'fontsize',8);hold on;
    xlabel('Frequency (Hz)','fontsize',8);
    ylabel({'Relative power'},'fontsize',8);
%     p=get(gca,'pos'); p(4)=p(4)+0.04; set(gca,'pos',p,'ylim',[-30 30]);
    
    ti_index=find(Events.II_time>w*60*s & Events.II_time<=w*60*(s+1));
    if isempty(ti_index), return; end
    fi=Fi(ti_index);        
    [y n] = gethistogram(fi,0.01);
    [freqi Pi]=FASPER(Events.II_time(ti_index),1./fi,4);
    h{1}=subplot(323);
    set(h{1},'fontsize',5);
    stem(n,1*y,'marker','o','markersize',2,'color','r'); ylim('auto'); 
    ylabel({'Probability'},'fontsize',8);hold on; 
%     plot(freqi,(Pi),'color',c2);hold on;
%     
%     p=get(gca,'pos');p(2)=p(2)-0.04;p(4)=p(4)+0.04;  set(gca,'pos',p);
    
    subplot(325)
    xlim([0 max([max(f);max(fi)])]/1);hold off; box off;
    subplot(323);   
    if isfield(Events,'UU_interval')
        tu_index=find(Events.UU_time>w*60*s & Events.UU_time<=w*60*(s+1));
        if isempty(tu_index),return; end
        fu=1./(Events.UU_interval(tu_index)); 
        if strcmp(label(1:3),'CMV'),dx=0.001;else dx=0.01;end
        [y n]=gethistogram(fu,dx);
        
        stem(n,1*y,'color',c1,'marker','^','markersize',2); 
%         ylim('auto'); %ylabel('histogram of Fu');
        legend('Breathing','Ventilator'); set(legend,'edgecolor',[1 1 1]);
    end    
    xlabel('Frequency (Hz)','fontsize',8);
    title({'InstFreq distribution of II&UU interval '},'fontsize',8);
    xlim([0 max([max(f);max(fi)])]/1);hold off; box off;

    % Plot Events
%     subplot(4,2,[1 2]),plot(tr(indexrr)/60,60./RR,'-','color','c'); axis 'tight'; hold on;
%     set(gca,'fontsize',5);
%     plot(tr(indexrr)/60,60./RR,'.','color','b'); axis 'tight';
%     ylabel('InstHR','fontsize',8);hold on;
%     xlabel('time (min)','fontsize',8);
%     xlim([tr(indexrr(1))/60 tr(indexrr(1))/60+1]);

%     temp=60*range(1./RR)/2.5;
%     hi=stem(Events.I_time(ti_index)/60,(60/max(RR)+temp)*ones(size(ti_index)),'r','marker','none','linewidth',1);
%     set(hi,'basevalue',60/max(RR));
%     if isfield(Events,'UU_interval')
%         hu=stem(Events.U_time(tu_index)/60,(60/min(RR)-temp)*ones(size(tu_index)),'color',c1,'marker','none','linewidth',1);
%         set(hu,'basevalue',60/min(RR));
%     end
%     
%     hr = stem(tr(indexrr)/60,(60/max(RR)+temp/5)*ones(size(indexrr)),'color',c2,'marker','none','linewidth',0.5);
%     set(hr,'basevalue',60/max(RR));
%     hold off;
%     ylim([60/max(RR) 60/min(RR)]);
%     title(figname,'fontsize',8);
end
