function [Indexall bsl rec xtickvec xlabelvec indices]= loadindex(filename,couptype,pairname,ismultivent,normalize)
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


Indexbsl=[];Indexvent=[];Indexrec=[];
for i=1:length(fileinfo.subjectlist)
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline',[filename '.mat'])    
    [Indexbsl name medianbsl]=getindex(Indexbsl,filepath,couptype,pairname,fileinfo,numGroup,normalize,[],i,0);
    Indexventi=[];
    for j=1:fileinfo.numvent(i)
        ventpart=['v' char(num2str(j))];        
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'ventilation',[filename ventpart '.mat']);
        [Indexventi name]=getindex(Indexventi,filepath,couptype,pairname,fileinfo,numGroup,normalize,medianbsl,i,j);
    end   
    Indexvent=[Indexvent;Indexventi];
    filepath = fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery',[filename '.mat']);        
    [Indexrec]=getindex(Indexrec,filepath,couptype,pairname,fileinfo,numGroup,normalize,medianbsl,i,j+1);
end
if ~ismultivent
    Indexbsl(:,end-1)=0;
    Indexvent(:,end-1)=1;
    Indexrec(:,end-1)=2;
end
indices = cell(length(name),1);
for i=1:length(name)
    indices{i}=[char(name(i)) '.' pairname];
end

nv = [0;3;5;7;9];
% nv =[0 1 2 3 4 5 6 7 8 9];
subIndexvent = cell(length(nv)-1,1);
for i=1:length(nv)-1
    subIndexvent{i} = Indexvent(Indexvent(:,end-1)>nv(i) & Indexvent(:,end-1)<=nv(i+1),:);
end

bsl = 0; xtickvec = bsl; xlabelvec={'bsl'};
v = zeros(length(nv)-1,1);
if ismultivent
    dv=0.75;     
    v(1) = bsl+1; xtickvec = [xtickvec v(1)]; %xlabelvec=[xlabelvec,['v' num2str(1)]]; 
    xlabelvec=[xlabelvec,'v1-3']; 
    text2={'4-5','6-7','8-9'};
    length(v)
    
    for i=2:length(v)
        i
        v(i) = v(i-1)+dv;
        xtickvec = [xtickvec v(i)]; xlabelvec=[xlabelvec,['v' text2{i-1}]];
    end
    rec=v(end)+1;  
    xtickvec = [xtickvec rec]; xlabelvec=[xlabelvec,'rec'];
else
    v=1; rec = 2;
    xtickvec=[bsl v rec]; xlabelvec={'bsl','vent','rec'};
end
Indexbsl=[Indexbsl bsl*ones(size(Indexbsl,1),1)];
for i=1:length(v)
    subIndexvent{i} = [subIndexvent{i} v(i)*ones(size(subIndexvent{i},1),1)];
end
Indexrec=[Indexrec rec*ones(size(Indexrec,1),1)];

Indexall = Indexbsl;
for i=1:length(v)
    Indexall=[Indexall;subIndexvent{i}];
end
Indexall=[Indexall;Indexrec]; 
ncol = size(Indexall,2);
Indexall = sortrows(Indexall,[ncol-4 ncol-2 ncol-3]);
% [index animalnumber epochnumber ventnumber bvv/cmv group(bsl/vent/rec] 

