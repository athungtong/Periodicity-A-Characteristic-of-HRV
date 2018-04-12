function [Indexall bsl rec xtickvec xlabelvec indices]= loadPoincareSDindex(filename,pairname,ismultivent,normalize,nump,Maxtau)
%% This function loads Index such as CVCIndex, RVCIndex, CRCIndex, and
% HRVIndex from the file "AllIndex.mat".
% [Index #Rat #Vent #Group #Period]
if nargin <3
    ismultivent=0;
    normalize=0;
end
eventpath = fullfile('CVC Events');
% load(fullfile(eventpath,'fileinfo2.mat'));
% numGroup=[zeros(6,1);ones(6,1)]; %1 for BVV 0 for CMV

load(fullfile(eventpath,'fileinfo4.mat'));
numGroup=[zeros(7,1);ones(8,1)]; %1 for BVV 0 for CMV
% numGroup=[zeros(1,1);ones(1,1)]; %1 for BVV 0 for CMV

w=1;
Indexbsl=cell(nump,Maxtau) ;Indexvent=cell(nump,Maxtau);Indexrec=cell(nump,Maxtau);
for i=1:length(fileinfo.subjectlist)
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline',[filename '.mat']);    
    [Indexbsl name medianbsl]=getPoincareSDindex(Indexbsl,pairname,filepath,fileinfo,numGroup,normalize,cell(nump,Maxtau),i,0,nump,Maxtau);
    Indexventi=cell(nump,Maxtau);
    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];        
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',[filename ventpart '.mat']);
        [Indexventi name]=getPoincareSDindex(Indexventi,pairname,filepath,fileinfo,numGroup,normalize,medianbsl,i,j,nump,Maxtau);
    end
    if ~isempty(Indexventi)
        for p=1:nump
            for t=1:Maxtau                   
%                 if ismultivent
%                     Indexventi{p,t}(:,end-1)=(w:w:size(Indexventi{p,t},1)*w)';
%                 end
                Indexvent{p,t}=[Indexvent{p,t};Indexventi{p,t}];
            end
        end
    end  
    
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery',[filename '.mat']);        
    [Indexrec]=getPoincareSDindex(Indexrec,pairname,filepath,fileinfo,numGroup,normalize,medianbsl,i,j,nump,Maxtau);
end

indices = cell(length(name),1);
for i=1:length(name)
    indices{i}=char(name(i));
end

bsl = 0;
if ismultivent
    dv=0.75;    
    v1=bsl+1; v2=v1+dv; v3 = v2+dv; v4 = v3+dv; 
    rec = v4+1;
    xtickvec=[bsl v1 v2 v3 v4 rec]; xlabelvec={'bsl','v1','v2','v3','v4','rec'};
else
    dv=0;
    v1=bsl+1; v2=v1+dv; v3 = v2+dv; v4 = v3+dv; 
    rec = v4+1;
    xtickvec=[bsl v1 rec]; xlabelvec={'bsl','vent','rec'};
end

for p=1:nump
    for t=1:Maxtau                 
%         nv1 = round(1/4*max(Indexvent{p,t}(:,end-1)));
%         nv2 = round(2/4*max(Indexvent{p,t}(:,end-1)));
%         nv3 = round(3/4*max(Indexvent{p,t}(:,end-1)));
%         nv4 = max(Indexvent{p,t}(:,end-1));
        nv1 = 3;    nv2=5;  nv3=7;  nv4=13;
        Indexvent1=Indexvent{p,t}(Indexvent{p,t}(:,end-1)>0 & Indexvent{p,t}(:,end-1)<=nv1,:);
        Indexvent2=Indexvent{p,t}(Indexvent{p,t}(:,end-1)>nv1 & Indexvent{p,t}(:,end-1)<=nv2,:);
        Indexvent3=Indexvent{p,t}(Indexvent{p,t}(:,end-1)>nv2 & Indexvent{p,t}(:,end-1)<=nv3,:);
        Indexvent4=Indexvent{p,t}(Indexvent{p,t}(:,end-1)>nv3 & Indexvent{p,t}(:,end-1)<=nv4,:);

        Indexbsl{p,t}=[Indexbsl{p,t} bsl*ones(size(Indexbsl{p,t},1),1)];
        Indexvent1=[Indexvent1 v1*ones(size(Indexvent1,1),1)];
        Indexvent2=[Indexvent2 v2*ones(size(Indexvent2,1),1)];
        Indexvent3=[Indexvent3 v3*ones(size(Indexvent3,1),1)];
        Indexvent4=[Indexvent4 v4*ones(size(Indexvent4,1),1)];
        Indexrec{p,t}=[Indexrec{p,t} rec*ones(size(Indexrec{p,t},1),1)];

        Indexall{p,t}=[Indexbsl{p,t};Indexvent1;Indexvent2;Indexvent3;Indexvent4;Indexrec{p,t}];
    end
end

end

function [Indexout name medianbsl]=getPoincareSDindex(Indexout,pairname,filepath,fileinfo,numGroup,normalize,medianbsl,i,j,nump,Maxtau)
%%
name{1} = 'sd2';name{2}='sd1';

if ~isfilefound(filepath),return;end
load(filepath);
if ~isfield(AllIndex,'TPVA'),return;end
if ~isfield(AllIndex.TPVA,pairname),return;end
for p=1:nump
    for t=1:Maxtau 
        index=[];
        index(:,1) = AllIndex.TPVA.(pairname).par(p).delay(t).sd2;   
        index(:,2) = AllIndex.TPVA.(pairname).par(p).delay(t).sd1; 
        if normalize && ~isempty(medianbsl{p,t})%normalize for vent and rec
            for n=1:length(name)
                index(:,n)=index(:,n)-medianbsl{p,t}(:,n); 
            end
        end
        if normalize && isempty(medianbsl{p,t})%normalize for bsl itself
            for n=1:length(name)
                medianbsl{p,t}(:,n) = median(index(:,n)); 
                index(:,n)=index(:,n)-medianbsl{p,t}(:,n);
            end    
        end
        Indexout{p,t}=[Indexout{p,t};[index fileinfo.subjectlist(i,:)*ones(size(index,1),1) j*ones(size(index,1),1) numGroup(i,:)*ones(size(index,1),1)]];           
    end
end
end
