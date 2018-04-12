function [PX PY PZ]=TimefrequencyPLOT(x,y,z,tx,ty,tz,fx,fy,fz,isplot,fname,PX,PY,PZ,freq)
% compute cross correlation between Lomb Periodogram of RR interval and
% histogram of Inspiration and Vent
% x,y are string specify 'I','U','R'
% tx,ty,fx,fy time and freq of x and y
% c is the maximum cross correlation


w=2.5;

Seg=floor(max(tx(end),ty(end))/60/w); %number of section, round up
if Seg==0,Seg=1;end


%freq = (0.1:0.001:2);
for s=0:Seg-1
   indextx= tx>=w*s*60 & tx<w*(s+1)*60;
   indexty= ty>=w*s*60 & ty<w*(s+1)*60;
   indextz= tz>=w*s*60 & tz<w*(s+1)*60;
   if isempty(indextx) || isempty(indexty),continue;end


   [nx Px]=FASPER(tx(indextx),1./fx(indextx));

   Px(nx<0.1)=0;

   Py = histc(fy(indexty),nx); 
   Py=Py/length(fy(indexty)); 
   ny=nx;

   % normalize to have 0-1 range
   Px=Px/sum(Px); Py=Py/sum(Py);
   %Py(Py==0)=NaN;
   if isempty(Py),continue; end

   %Compute Linear interpolation for new nx, ny
    a = [1 0.2];
    b = [2 3];
     Px = smooth(Px);
     %plot(ny,Py); hold on;
     Py = smooth(Py);    
     %plot(ny,Py,'r'); hold off; input(':');
    
   Px = interp1(nx,Px,freq);
   Py = interp1(ny,Py,freq);
   %Px = (Px-min(Px))/(range(Px))*255;
   %Py = (Py-min(Py))/(range(Py))*255;
   PX=[PX Px'];
   PY=[PY Py'];

   if ~isempty(indextz)
       Pz = histc(fz(indextz),nx); Pz=Pz/length(fz(indextz)); nz=nx;

       Pz=Pz/sum(Pz);
       %Pz(Pz==0)=NaN;
       if isempty(Pz),continue; end
       Pz = smooth(Pz);
       
       Pz = interp1(nz,Pz,freq);
       
       %Pz = (Pz-min(Pz))/(range(Pz))*255;
       PZ=[PZ Pz'];
   end

end