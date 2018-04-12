
% 1-7 [MNN SDNN Flomb.LFP Flomb.HFP Flomb.LHR DFA Poincare.SD1 ... 
% 8-11                Poincare.SD2 Poincare.SDRatio Poincare.UP Poincare.DOWN ...
% 12-14               DelayedPoincare.SD1 DelayedPoincare.SD2 DelayedPoincare.SDRatio...
% 15-20                DelayedPoincare.UP DelayedPoincare.DOWN CE DelayedCE (s+1)*w w]; 
   
                
function tempMergeIndex(indexname,pairname)
%%
drive = fullfile('D:');
eventpath = fullfile(drive,'Dropbox','Research','Matt CVC project','CVC Events');
load(fullfile(eventpath,'fileinfo.mat'));


for i=17:length(fileinfo.subjectlist)
    
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','AllIndex1.mat');  
    if isfilefound(filepath)
        temp1=load(filepath);
        filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','AllIndex2.mat');    
        temp2=load(filepath);
        
        name=fieldnames(temp1.AllIndex.(indexname).(pairname));
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'baseline','AllIndex.mat');
        if isfilefound(filepath)
            load(filepath);
            if isfield(AllIndex,indexname)
                if isfield(AllIndex.(indexname),pairname)
                    AllIndex.(indexname) = rmfield(AllIndex.(indexname),pairname);
                end
            end            
        else
            AllIndex=struct;
        end
%         temp2.AllIndex.HRV.RR
        for j=1:length(name)
            Index1=temp1.AllIndex.(indexname).(pairname).(char(name(j)));
            Index2=temp2.AllIndex.(indexname).(pairname).(char(name(j)));
            AllIndex.(indexname).(pairname).(char(name(j)))=[Index1;Index2];%(Index1+Index2)/2;%
        end        
        save(filepath,'AllIndex');        
    end
    
    %%%%%%%%%%
    filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','AllIndex1.mat');   
    if isfilefound(filepath)
        temp1=load(filepath);
        filepath= fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','AllIndex2.mat');    
        temp2=load(filepath);
        
        name=fieldnames(temp1.AllIndex.(indexname).(pairname));
        filepath=fullfile(eventpath,fileinfo.Group(i,:),fileinfo.Subject(i,:),'recovery','AllIndex.mat');
        if isfilefound(filepath)
            load(filepath);
            if isfield(AllIndex,indexname)
                if isfield(AllIndex.(indexname),pairname)
                    AllIndex.(indexname) = rmfield(AllIndex.(indexname),pairname);
                end
            end            
        else
            AllIndex=struct;
        end

        for j=1:length(name)
            Index1=temp1.AllIndex.(indexname).(pairname).(char(name(j)));
            Index2=temp2.AllIndex.(indexname).(pairname).(char(name(j)));
            AllIndex.(indexname).(pairname).(char(name(j)))=[Index1;Index2];%(Index1+Index2)/2;%
        end        
        save(filepath,'AllIndex');        
    end
end


