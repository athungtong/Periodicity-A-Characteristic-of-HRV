function [Indexout name medianbsl]=getindex(Indexout,filepath,couptype,pairname,fileinfo,numGroup,normalize,medianbsl,i,j)
%%
name=[];
if ~isfilefound(filepath),return;end
load(filepath);
if ~isfield(AllIndex,couptype),return;end
if ~isfield(AllIndex.(couptype),pairname),return;end

temp=getfield(AllIndex.(couptype),pairname);
name = fieldnames(temp);    index=[];
lnindex = min(structfun(@length,AllIndex.(couptype).(pairname)));
index = zeros(lnindex,length(name));
for n=1:length(name)
    nextindex = getfield(temp,char(name(n)));
    index(1:lnindex,n)= nextindex(1:lnindex);
end

if normalize && ~isempty(medianbsl)%normalize for vent and rec
    for n=1:length(name)
        index(:,n)=index(:,n)-medianbsl(:,n); 
    end
end

if normalize && isempty(medianbsl)%normalie for bsl itself
    for n=1:length(name)
        medianbsl(:,n) = median(index(:,n)); 
        index(:,n)=index(:,n)-medianbsl(:,n);
    end    
end

% if normalize,index=index-ones(size(index,1),1)*median(index);end
%need to change j to reflect real vent number
if fileinfo.subjectlist(i) == 27
    if j>1
        j=j-1;
    end
end
if fileinfo.subjectlist(i) == 25
    if j == 2 || j==3 || j==4,j=2;end
    if j==5 || j==6, j=3;end
    if j==9 || j==10,j=6;end
    if j>=11,j=j-4;end
end
% size(index)
    Indexout=[Indexout;[index fileinfo.subjectlist(i,:)*ones(size(index,1),1) (1:size(index,1))' j*ones(size(index,1),1) numGroup(i,:)*ones(size(index,1),1)]];  %put epoch num in         
% [index animalnumber epochnumber ventnumber bvv/cmv] 