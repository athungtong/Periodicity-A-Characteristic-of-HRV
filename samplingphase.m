function [t phx2 phy2]=samplingphase(tx,ty,phx,phy,tau)

t1=max(tx(1),ty(1));
tn = min(tx(end),ty(end));

t=(t1:tau:tn)';
phx2=interp1(tx,phx,t);
phy2=interp1(ty,phy,t);

phx2(isnan(phx2))=[];
phy2(isnan(phy2))=[];

% [txy phx2 phy2] = phaserelative(tx,ty,0);
% t=(txy(1):tau:txy(end))';    
% phx2=interp1(txy,phx2,t);
% phy2=interp1(txy,phy2,t);
% phx2(isnan(phx2))=[];
% phy2(isnan(phy2))=[];
