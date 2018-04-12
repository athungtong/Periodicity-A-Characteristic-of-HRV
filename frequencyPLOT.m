function c=frequencyPLOT(x,y,z,tx,ty,tz,fx,fy,fz,h1,h2,h3,fname)
% compute cross correlation between Lomb Periodogram of RR interval and
% histogram of Inspiration and Vent
% x,y are string specify 'I','U','R'
% tx,ty,fx,fy time and freq of x and y
% c is the maximum cross correlation


w=2.5;

Seg=floor(max(tx(end),ty(end))/60/w); %number of section, round up
if Seg==0,Seg=1;end
c=-10*ones(Seg,1);
for s=0:Seg-1
    s=2;
   indextx= tx>=w*s*60 & tx<w*(s+1)*60;
   indexty= ty>=w*s*60 & ty<w*(s+1)*60;
   indextz= tz>=w*s*60 & tz<w*(s+1)*60;
   if isempty(indextx) || isempty(indexty),continue;end

   xisR=0; yisR=0; 
       [nx Px]=FASPER(tx(indextx),1./fx(indextx));
       Px(nx<0.1)=0;
       Py = histc(fy(indexty),nx); Py=Py/length(fy(indexty)); ny=nx;
           
       % normalize to have 0-1 range
       Px=Px/sum(Px); Py=Py/sum(Py);

       Py(Py==0)=NaN;
       if isempty(Py),continue; end
       
       axes(h1)
       stem(ny(1:1:end),Py(1:1:end),'color',[0,99,255]/255,'markersize',2.0,'ShowBaseLine','off');hold on;
       box off;
       
       axes(h2), plot(nx,Px,'color','r','markersize',2.5);
       box off;
       %ylabel('Power of RR-Interval'); hold off;
      
       if ~isempty(indextz)
           Pz = histc(fz(indextz),nx); Pz=Pz/length(fz(indextz)); nz=nx;
           Pz=Pz/sum(Pz);
           Pz(Pz==0)=NaN;
           axes(h1), stem(nz(1:1:end),Pz(1:1:end),'^','color',[173,118,0]/225,'markersize',2.0,'ShowBaseLine','off');hold off;
       end
       
       axes(h3)
       
       gamma=acf(fx(indextx),100);
       plot(gamma(2:end),'color',[0,102,51]/255); xlim([0 100]); ylim([-1 1]);
       box off;
       %tight_subplot(2,1,[.01 .03]);
       %print( fig, '-dpdf', [fname2 '.pdf'] );
       %input('continue?');
       
       break; 
   
      maxlag = round(0.05*length(Px));
      c(s+1)=max(abs(xcorr(center(Px),center(Py),maxlag,'unbiased')));

end
return;
c(c==-10)=[];