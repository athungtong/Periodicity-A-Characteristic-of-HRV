function pvalue = mainmanova(x,y,h,pos,fid,xlabel)
%%
c1=[0 100 255]/255;
c2=[200 0 0]/255;



numratx=unique(x(:,2));
numepochx=zeros(size(numratx));
for i=1:length(numratx)
    numepochx(i)=length(x(x(:,2)==numratx(i)));
end
[numratx numepochx];

numraty=unique(y(:,2));
numepochy=zeros(size(numraty));
for i=1:length(numraty)
    numepochy(i)=length(y(y(:,2)==numraty(i)));
end
[numraty numepochy];
minepoch = min([numepochx;numepochy]);

% method 1, use the smallest number of epoch
tempx=[];
for i=1:length(numratx)
    temp=x(x(:,2)==numratx(i),:);
    tempx = [tempx;(temp(randperm(numepochx(i),minepoch),1))'];
end


tempy=[];
for i=1:length(numraty)
    temp=y(y(:,2)==numraty(i),:);
    tempy = [tempy;(temp(randperm(numepochy(i),minepoch),1))'];
end

[p, table] = anova_rm({tempx tempy}, 'off');


xdat=reshape(tempx,size(tempx,1)*size(tempx,2),1);
temp=prctile(xdat,[25 50 75]);
p25x = temp(1);
p50x = temp(2);
p75x = temp(3);


ydat=reshape(tempy,size(tempy,1)*size(tempy,2),1);
temp=prctile(ydat,[25 50 75]);
p25y = temp(1);
p50y = temp(2);
p75y = temp(3);

pvalue = p(2);

%  data =[zeros(length(xdat),1) xdat;ones(length(ydat),1) ydat];
%  xlswrite('MIX.xls',data);


boxplot(xdat,'position',pos,'color',c1,'notch','on','widths',0.2,'symbol','.');hold on;
set(h,'XTickLabel',{''})
set(h,'XTickLabel',{''});
h2 = findobj(gca,'Tag','Box');
for j=1:2:length(h2)
patch(get(h2(j),'XData'),get(h2(j),'YData'),[128,229,255]/255,'FaceAlpha',.5);
end
%lines = findobj(gca, 'type', 'line', 'Tag', 'Median');
%set(lines,'linewidth',1, 'Color', [128,229,255]/255); 
%  
 
boxplot(ydat,'position',pos+.22,'color',c2,'notch','on','widths',0.2,'symbol','.');hold on;
xlim([pos(1)-0.5 pos(end)+0.5]);
set(h,'XTickLabel',{''})



fprintf(fid,'%s\r\n',['    During: ' xlabel]);
fprintf(fid,'\t%s\r\n',['BVV animal #: ',num2str(numratx')]);
fprintf(fid,'\t%s\r\n',['CMV animal #: ',num2str(numraty')]);
fprintf(fid,'\t%s\r\n',[num2str(length(numratx)*minepoch) ' BVV,' num2str(length(numraty)*minepoch) ' CMV, (' num2str(minepoch) ' each)' ]);
fprintf(fid,'\t%s\r\n',['BVV stat: 25,50,75 percentiles: ', num2str([p25x p50x p75x])]);
fprintf(fid,'\t%s\r\n',['CMV stat: 25,50,75 percentiles: ', num2str([p25y p50y p75y])]);
fprintf(fid,'\t%s\r\n\r\n',['p-value: ',num2str(pvalue)]);


% manova1([tempx;tempy(1:7,:)],[groupx;groupy(1:end-1)])
% axis(h);
% bar(pos,p50x,'barwidth',0.3,'facecolor',[1 1 1]);hold on;
% bar(pos+0.32,p50y,'barwidth',0.3,'facecolor',[0.494 0.494 0.494]);
% errorbar(pos,p50x,p50x-p25x,p75x-p50x,'color','k');
% errorbar(pos+0.32,p50y,p50y-p25y,p75y-p50y,'color','k');

