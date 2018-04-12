function mainCRCstat(filename,couptype,pairname,ismultivent,normalize,numindex,numplot)
%%

[CVCall bsl rec xtickvec xlabelvec indices]=loadindex(filename,couptype,pairname,ismultivent,normalize);
% Transformation
% CVCall(CVCall(:,3)==Inf,3)=1000; 
% CVCall(CVCall(:,4)==Inf,4)=1000; 
% CVCall(:,1:4)=log(CVCall(:,1:4));

CVCBVV=CVCall(CVCall(:,end-1)==0,:);
CVCCMV=CVCall(CVCall(:,end-1)==1,:);
%%
figure(1);
numpair = 4;
for i=1:numindex
    subplot(numindex,numpair,numplot+numpair*(i-1))
    boxplot(CVCBVV(:,i),CVCBVV(:,end),'position',CVCBVV(:,end),'color','b','notch','on','widths',0.2,'symbol','b+'); hold on;
    set(gca,'XTickLabel',{' '})

    boxplot(CVCCMV(:,i),CVCCMV(:,end),'position',CVCCMV(:,end)+.22,'color','r','widths',0.2,'notch','on','symbol','r+'); hold off;
    set(gca,'XTick',xtickvec+.11,'XTickLabel',xlabelvec)
    ylim 'auto'; title(indices{i});hold off; 
    xlim([bsl-0.25 rec+0.5]);
    grid on
    
    x=[CVCBVV(:,numindex) CVCBVV(:,end)]; x(x(:,1)<0,:)=[]; %x(:,1)=sqrt(x(:,1));
    x1 = x(x(:,2)==0); x2 = x(x(:,2)==1); x3 = x(x(:,2)==2);
    [mean(x1) mean(x2) mean(x3);std(x1) std(x2) std(x3)]

    y=[CVCCMV(:,numindex) CVCCMV(:,end)]; y(y(:,1)<0,:)=[]; %y(:,1)=sqrt(y(:,1));
    y1 = y(y(:,2)==0); y2 = y(y(:,2)==1); y3 = y(y(:,2)==2);
    [mean(y1) mean(y2) mean(y3);std(y1) std(y2) std(y3)]
%     [h p]=ranksum(x1,y1,0.05); [h p]
    [h p]=ranksum(x2,y2,0.05); [h p]
%     [h p]=ranksum(x3,y3,0.05); [h p]
    
%     [h p]=ttest2(sqrt(x1),sqrt(y1),0.05); [h p]
    [h p]=ttest2(sqrt(x2),sqrt(y2),0.05); [h p]
%     [h p]=ttest2(sqrt(x3),sqrt(y3),0.05); [h p]
% save('BVVI.txt','x','-ascii');
% save('CMVI.txt','y','-ascii');

end    

