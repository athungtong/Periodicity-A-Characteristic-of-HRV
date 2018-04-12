function tend=mainFreqmod(filepath,fOutName,technique,w,tau,couptype,tend,h1,h2,h3,fname)


tend=0;

if ~isfilefound(filepath),return;end
subjectpath = fileparts(filepath);
load(filepath);
% Check breath. For all coupling type, we evaluate for at least 200
% breath
% numbreath=200;
% if length(Events.I_time) < numbreath
%     display(['Number of breath is less than ' num2str(numbreath)]);
%     return;
% end


x='R';y='I';z='U'; 

if ~isfield(Events,[x x '_time']),return;end 
if ~isfield(Events,[y y '_time']),return;end

if ~isfield(Events,[z z '_time'])
    tz=-1;
    fz=-1;
else
   tz=getfield(Events,[z z '_time']);
   fz=1./getfield(Events,[z z '_interval']);
end

tx=getfield(Events,[x x '_time']);
% save([fname '.txt'],'tx','-ascii');
% return
ty=getfield(Events,[y y '_time']);

fx=1./getfield(Events,[x x '_interval']);
fy=1./getfield(Events,[y y '_interval']);

c=frequencyPLOT(x,y,z,tx,ty,tz,fx,fy,fz,h1,h2,h3,fname);

return;

if isfilefound(fullfile(subjectpath,fOutName))
   load(fullfile(subjectpath,fOutName)); %This will load AllIndex
else
    AllIndex=struct;
end

if isfield(AllIndex,'FMOD')
    if isfield(AllIndex.FMOD,([x y]))
        AllIndex.FMOD = rmfield(AllIndex.FMOD,([x y]));
    end
end
  
AllIndex.FMOD.([x y]).freqx=c;

if ~isplot
    save(fullfile(subjectpath,fOutName),'AllIndex');
end



