function [Ds I12s I21s MIs]=direcsurrogate_time(tx,ty,safr,tau,nbin,Ns,fignum,theta1,theta2,Dtheta1,Dtheta2,isplot)

%% Do surrogate test
if nargin < 12
    isplot=0;
end
Ds=zeros(Ns,1); I12s=zeros(Ns,1);   I21s=zeros(Ns,1);  MIs=zeros(Ns,1);
% phx = 2*pi*(0:length(tx)-1)';
% phy = 2*pi*(0:length(ty)-1)'; 
% [~, phx phy]=phaserelative(tx,ty,1);
% cindex=1-circularvar(phasediff(phx,phy,1,1));

% figure(5);stroboscope(tx,ty,1);
for k=1:Ns
    [txs tys]=phaseSurrogate2(tx,ty);
    pxs=(0:length(txs)-1)'*2*pi;
    pys=(0:length(tys)-1)'*2*pi;
    [ts, phxs phys]=samplingphase(txs,tys,pxs,pys,1/safr);
    [Ds(k) I12s(k) I21s(k) theta1s theta2s Dtheta1s Dtheta2s MIs(k)]=main_coupdirecinfo(phxs,phys,safr,tau,nbin);

%%     figure(2)
%     subplot(211);plot(txs(1:end-1),diff(txs),'r.-');hold on;
%     subplot(212);plot(tys(1:end-1),diff(tys),'r.-');hold on;
    
%     phxs = 2*pi*(0:length(txs)-1)';
%     phys = 2*pi*(0:length(tys)-1)';

%     figure(4);stroboscope(txs,tys,1);
%     input('con');
%     [ts';mod(phxs,2*pi)';mod(phys,2*pi)']
    %     tys=phaseSurrogate(ty,tx,phy,phx);
%     [~, phxs phy]=samplingphase(tx,tys,1/safr);

%     figure(3)
%     subplot(121);plot(ts(1:end),mod(phxs,2*pi),'r.-');hold on;
%     subplot(122);plot(ts(1:end),mod(phys,2*pi),'r.-');hold on;
%%
%     I12si=zeros(size(tau));      I21si=zeros(size(tau));
%     MIsi=I12si;
%     for i=1:length(tau)
%         [theta1s,theta2s,Dtheta1s,Dtheta2s]=getDtheta(phxs,phys,1,safr,tau(i));
%         [~, I12si(i) I21si(i) MIsi(i)] = coupdirecinfo(mod(theta1s,2*pi),mod(theta2s,2*pi),Dtheta1s,Dtheta2s,2,nbin); 
% %         I12si(i) = condmutualinfo(mod(theta1s,2*pi),Dtheta2s,mod(theta2s,2*pi),2,nbin); 
% %         I21si(i) = condmutualinfo(mod(theta2s,2*pi),Dtheta1s,mod(theta1s,2*pi),2,nbin);
%     end
%     I12s(k) = mean(I12si);
%     I21s(k) = mean(I21si);
%     MIs(k) = mean(MIsi);
% %     Ds(k) = mean( (I12si-I21si)./(I12si+I21si) );
%     Ds(k) = (I12s(k)-I21s(k))/(I12s(k)+I21s(k));
    
%     figure(1);
%     subplot(211);plot(Dtheta1s,'r.-');hold on;
%     subplot(212);plot(Dtheta2s,'r.-');hold on;

    if isplot
        figure(fignum)
        subplot(241); [f,xi] = ksdensity(mod(theta1s,2*pi)); plot(xi,f,'r'); ylabel('phi xs'); hold on;
        subplot(242); [f,xi] = ksdensity(mod(theta2s,2*pi)); plot(xi,f,'r'); ylabel('phi ys');hold on;
        subplot(243); [f,xi] = ksdensity(Dtheta1s); plot(xi,f,'r'); ylabel('Dphi xs');hold on;
        subplot(244); [f,xi] = ksdensity(Dtheta2s); plot(xi,f,'r' ); ylabel('Dphi ys');hold on;
    end
end
if isplot
    figure(fignum)
    subplot(241); [f,xi] = ksdensity(mod(theta1,2*pi)); plot(xi,f); ylabel('phi x'); hold off;
    subplot(242); [f,xi] = ksdensity(mod(theta2,2*pi)); plot(xi,f); ylabel('phi y');hold off;
    subplot(243); [f,xi] = ksdensity(Dtheta1); plot(xi,f); ylabel('Dphi x');hold off;
    subplot(244); [f,xi] = ksdensity(Dtheta2); plot(xi,f); ylabel('Dphi y');hold off;

    % subplot(245);X=getspectrum(mod(theta1s,2*pi)); plot(log(X(1:round(length(X)/2))),'r');hold on;
    % subplot(245);X=getspectrum(mod(theta1,2*pi)); plot(log(X(1:round(length(X)/2)))); hold off;
    % 
    % subplot(246);X=getspectrum(mod(theta2s,2*pi)); plot(log(X(1:round(length(X)/2))),'r');hold on;
    % subplot(246);X=getspectrum(mod(theta2,2*pi)); plot(log(X(1:round(length(X)/2)))); hold off;
    subplot(245);plot(mod(theta1(100:200),2*pi));hold on; plot(mod(theta1s(100:200),2*pi),'r'); hold off;
    subplot(246);plot(mod(theta2(100:200),2*pi));hold on; plot(mod(theta2s(100:200),2*pi),'r'); hold off;

    % subplot(247);X=getspectrum(Dtheta1s); plot(log(X(1:round(length(X)/2))),'r');hold on;
    % subplot(247);X=getspectrum(Dtheta1); plot(log(X(1:round(length(X)/2)))); hold off;
    % 
    % subplot(248);X=getspectrum(Dtheta2s); plot(log(X(1:round(length(X)/2))),'r');hold on;
    % subplot(248);X=getspectrum(Dtheta2); plot(log(X(1:round(length(X)/2)))); hold off;
    subplot(247);plot(Dtheta1);hold on; plot(Dtheta1s,'r'); hold off;
    subplot(248);plot(Dtheta2);hold on; plot(Dtheta2s,'r'); hold off;
    % figure(1);
    % subplot(211); plot(Dtheta1,'.-');hold off;
    % subplot(212); plot(Dtheta2,'.-');hold off;
    % 
    % figure(2)
    % subplot(211);plot(tx(1:end-1),diff(tx),'.-');hold off;
    % subplot(212);plot(ty(1:end-1),diff(ty),'.-');hold off;
    % 
    % figure(3)
    % subplot(121);plot(t(1:end),mod(phx2,2*pi),'.-');hold off;
    % subplot(122);plot(t(1:end),mod(phy2,2*pi),'.-');hold off;

    % figure(6);stroboscope(txs,tys,1);
end





















%%
% for k=1:Ns
%     txs=phaseSurrogate(tx,ty,phx,phy);
%     [~, phxs phy]=samplingphase(txs,ty,1/safr);
% %     tys=phaseSurrogate(ty,tx,phy,phx);
% %     [~, phxs phy]=samplingphase(tx,tys,1/safr);
% 
%     I12si=zeros(size(tau));
%     for i=1:length(tau)
%         [theta1s,theta2,~,Dtheta2]=getDtheta(phxs,phy,1,safr,tau(i));
%         I12si(i) = condmutualinfo(mod(theta1s,2*pi),Dtheta2,mod(theta2,2*pi),2,nbin); 
%     end
%     I12s(k) = mean(I12si);
% 
%     tys=phaseSurrogate(ty,tx,phy,phx);
%     [~, phx phys]=samplingphase(tx,tys,1/safr);
% %     txs=phaseSurrogate(tx,ty,phx,phy);
% %     [~, phx phys]=samplingphase(txs,ty,1/safr);
%     I21si=zeros(size(tau));
%     for i=1:length(tau)
%         [theta1,theta2s,Dtheta1]=getDtheta(phx,phys,1,safr,tau(i));
%         I21si(i) = condmutualinfo(mod(theta2s,2*pi),Dtheta1,mod(theta1,2*pi),2,nbin); 
%     end
%     I21s(k) = mean(I21si);
% end
