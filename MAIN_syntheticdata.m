%clear all
Mu = 1; w1 = 40;

max_i=3; sb=1; n=1;
dw = (0:0.05:1)*50;

da = (1:0.5:11)/1;
for i = 1:20
    for j=1:length(dw);

        [XY time]=vanderpol1(Mu,w1,dw(j));
        figure(1);outfig=3;
 
        [f P] = FASPER(time,XY(:,1),1);
        f(P<0)=[];  P(P<0)=[];
        relatP = P/sum(P);
        maxPER(i,j) = max(relatP);

        bin=20;
        dx = range(f)/bin;
        edges = min(f):dx:max(f);   
        edges(end)=Inf;
        P2=zeros(bin,1);
        for b = 1:bin
            P2(b) = sum(relatP(f>=edges(b) & f<edges(b+1)));
        end
        edges(P2==0)=[];
        P2(P2==0)=[];
        tRSE(i,j) = 1+dot(P2,log(P2))/log(bin);
        
        [acy dy]= acf(XY(:,1),200); 
        dy=dy(202:end); acy=acy(202:end);
        maxACF(i,j) = max(abs(acy(26:50)));
        meanACF(i,j) = mean(abs(acy(26:50)));
        
        gammar = acy;

        [~, t] = time_delayed_MI(gammar,0:10);
        [m maxnorm] = fnn(gammar,'max',t,5,0);
         SE(i,j)= Entropy(maxnorm,center(gammar),m,t,.3,'max');   

        gammar = XY(:,1);
        [~, t] = time_delayed_MI(gammar,0:10);
        [m maxnorm] = fnn(gammar,'max',t,5,0);
         SEraw(i,j)= Entropy(maxnorm,center(gammar),m,t,.3,'max');   

        if j==5 || j==10 || j==15
            figure(1)
            switch j
                case{5},sub=1;
                case{10},sub=2;
                case{15},sub=3;
            end
            h1=subplot(3,3,sub); plot(time,XY(:,1),'-','markersize',4,'color',[200 0 0]/255);
            xlim([0 15]);  set(gca,'fontsize',6);
            ylim([-4 4]);  set(gca,'ytick',[]);
%             p=get(h1,'pos');p(3:4) = p(3:4)+0.035; set(h1,'pos',p); 
            box off
            xlabel('time (s)','fontsize',8,'interpreter','latex');
           if j == 5
                ylabel('$x$','fontsize',8,'interpreter','latex'); set(gca,'ytickmode','auto');
                ylabh=get(h1,'ylabel');
                set(ylabh,'unit','normalize');
                ylabpos = get(ylabh,'pos');
                ylabpos(1) = -0.3;
                set(ylabh,'pos',ylabpos);
           end    
            title(['$C= $' num2str(dw(j))],'fontsize',8,'interpreter','latex');
            tit=get(h1,'title'); set(tit,'unit','normalize');
            titpos=get(tit,'pos'); titpos(3)=-0.5;  set(tit,'pos',titpos);
            
            h2=subplot(3,3,3+sub);plot(f,(relatP),'color',[0 100 255]/255); 
            xlim([0.5 1.5]); set(gca,'fontsize',6); set(gca,'ytick',[]);box off
            xlabel('Frequency (Hz)','fontsize',8,'interpreter','latex');
            ylim([0 0.2]);
%             ylim([-15 0]);
%             p=get(h2,'pos'); p(3:4) = p(3:4)+0.035; set(h2,'pos',p)            
            if j == 5
                ylabel('Relative Power','fontsize',8,'interpreter','latex');
                set(gca,'ytickmode','auto');
                ylabh=get(h2,'ylabel');
                set(ylabh,'unit','normalize');
                ylabpos = get(ylabh,'pos');
                ylabpos(1) = -0.3;
                set(ylabh,'pos',ylabpos);

            end
            % ylim([0 0.5]); 
            h3=subplot(3,3,6+sub);plot(dy,acy,'-','markersize',4,'color',[0 100 0]/255);
            hold off; ylim([-1 1]); xlim([0 50]); set(gca,'fontsize',6);
             
            set(gca,'ytick',[]);box off
            xlabel('Delay (samples)','fontsize',8,'interpreter','latex');
%             p=get(h3,'pos'); p(3:4) = p(3:4)+0.035; set(h3,'pos',p)
            if j == 5
                ylabel('Autocorrelation','fontsize',8,'interpreter','latex');
                set(gca,'ytickmode','auto');
                ylabh=get(h3,'ylabel');
                set(ylabh,'unit','normalize');
                ylabpos = get(ylabh,'pos');
                ylabpos(1) = -0.3;
                set(ylabh,'pos',ylabpos);
            end
            % subplot(4,max_i,i+15),plot3(acy(1:end-2,1),acy(2:end-1,1),acy(3:end,1),'c-','linewidth',0.2);hold on;
            % subplot(4,max_i,i+15),plot3(acy(1:end-2,1),acy(2:end-1,1),acy(3:end,1),'b.','markersize',5);hold off;
            % subplot(4,max_i,sb+3*max_i),plot(acy(1:end-1,1),acy(2:end,1),'c-','linewidth',0.2);hold on;
            % subplot(4,max_i,sb+3*max_i),plot(acy(1:end-1,1),acy(2:end,1),'b.','markersize',5);hold off;
            % xlim([-1 1]); ylim([-1 1]);
        end

    end
    i
end

% fig1h=figure(1);
% set(fig1h,'paperposition',[0.25 2.5 3.417 2.8]);
% print(fig1h,'exVDP','-depsc');

%%
deg = zeros(5,1);
deg(1) = abs(degreeofmono(tRSE));
deg(2) = abs(degreeofmono(maxPER));
deg(3) = abs(degreeofmono(maxACF));
deg(4) = abs(degreeofmono(meanACF));
deg(5) = abs(degreeofmono(SE));
% deg(6) = degreeofmono(SEraw);

clear stattRSE statmaxPER statmaxACF statSE statSEraw;
for i=1:size(tRSE,2)
    stattRSE(i,:) = prctile(tRSE(:,i),[25 50 75]);
    statmaxPER(i,:) = prctile(maxPER(:,i),[25 50 75]);
    statmaxACF(i,:) = prctile(maxACF(:,i),[25 50 75]);
    statmeanACF(i,:) = prctile(meanACF(:,i),[25 50 75]);
    statSE(i,:) = prctile(SE(:,i)/2.5,[25 50 75]);
%     statSEraw(i,:) = prctile(SEraw(:,i),[25 50 75]);
end
%%

figure(3);
errorbar(dw,stattRSE(:,2),stattRSE(:,2)-stattRSE(:,1),stattRSE(:,3)-stattRSE(:,2),'marker','.','markersize',8,'color','k','linestyle','-');hold on;
errorbar(dw,statmaxPER(:,2),statmaxPER(:,2)-statmaxPER(:,1),statmaxPER(:,3)-statmaxPER(:,2),'marker','^','markersize',3,'color',[200 120 0]/255,'linestyle',':');hold on;
errorbar(dw,statmaxACF(:,2),statmaxACF(:,2)-statmaxACF(:,1),statmaxACF(:,3)-statmaxACF(:,2),'marker','s','markersize',3,'color',[0 100 255]/255,'linestyle',':');hold on;
errorbar(dw,statmeanACF(:,2),statmeanACF(:,2)-statmeanACF(:,1),statmeanACF(:,3)-statmeanACF(:,2),'marker','*','markersize',4,'color',[200 0 0]/255,'linestyle',':');hold on;


errorbar(dw,statSE(:,2),statSE(:,2)-statSE(:,1),statSE(:,3)-statSE(:,2),'marker','o','markersize',3,'color',[0 100 0]/255,'linestyle',':');hold on;
% errorbar(dw,statSEraw(:,2),statSEraw(:,2)-statSEraw(:,1),statSEraw(:,3)-statSEraw(:,2),'marker','x','markersize',4,'color',[150 0 150]/255,'linestyle',':');hold on;
xlim([-1 51]);
hold off;
set(gca,'fontsize',6);
xlabel('$C$ (Gaussian noise)','fontsize',8,'interpreter','latex');
ylabel('Proposed Measures of Periodicity','fontsize',8,'interpreter','latex');
legh=legend('tRSE','maxPER','maxACF','meanACF','SampEn');
set(legh,'Box','off','color','none','Position',[0.42,0.85,0.25,0.1],'fontsize',8,'interpreter','latex')
ylim([0 1.0]);
box off

text(30,1,'Degree of','fontsize',8,'interpreter','latex');
text(30,0.95,'Monotonicity','fontsize',8,'interpreter','latex');
xlabelvec =[{'tRSE'},{'maxPER'},{'maxACF'},{'meanACF'},{'SampEn'}];
pos = 0.85:-0.04:0;
for i=1:length(xlabelvec)
    text(30,pos(i),[num2str(floor(deg(i)*1000)/1000)],'fontsize',8,'interpreter','latex');
end

% fig1h=figure(1);
% set(fig1h,'paperposition',[0.25 2.5 4.3 2.8]);
% print(fig1h,'degVDP','-depsc');

