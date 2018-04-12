function [D I12 I21 MI]=coupdirecinfo(x,y,dx,dy,margin,par)

% % easier command but need more time
% %     I12 = condmutualinfo(x,dx,y);
% %     I21 = condmutualinfo(y,dy,x);
% %     D = (I12-I21)/(I12+I21);
% %
% % refference:
% % Bahraminasab, Alireza,Ghasemi, Fatemeh,Stefanovska, Aneta,McClintock, Peter V E, Kantz, Holger
% % Direction of coupling from phases of interacting oscillators: A permutation information approach.
% % Physical Review Letters  v. 100 issue 8, 2008
% 
% I12 = mutual_info(x,dy);
% I21 = mutual_info(y,dx);
% D=(I12-I21)/(I12+I21);
% return;
%% General marginal probability method
if margin == 1
    N=size(x,1);
    if nargin == 5
        par = round(sqrt(N));
    end
    nbin=par;
    [~, ~, atxbin]= get_margin_prob(x,nbin);
    [~, ~, atybin]= get_margin_prob(y,nbin);
    [~, ~, atdxbin]= get_margin_prob(dx,nbin);
    [~, ~, atdybin]= get_margin_prob(dy,nbin);
end

%% Marginal equiquantization method
if margin == 2
    if nargin == 5
        par = 16;
    end
    N=size(x,1); 
    [~, ~, atxbin]= get_margin_eq_prob(x,par);
    [~, ~, atybin]= get_margin_eq_prob(y,par);
    [~, ~, atdxbin]= get_margin_eq_prob(dx,par);
    [~, nbin, atdybin]= get_margin_eq_prob(dy,par);
end
% figure(8)
% [f,xi] = ksdensity(x); 
% subplot(2,2,1),plot(xi,f);
% [f,xi] = ksdensity(y); 
% subplot(2,2,2),plot(xi,f);
% [f,xi] = ksdensity(dx); 
% subplot(2,2,3),plot(xi,f);
% [f,xi] = ksdensity(dy); 
% subplot(2,2,4),plot(xi,f);




%% Permutation probability method
if margin == 3
    if nargin == 5
        par = 2;
    end
    N=size(x,1)-par+1;
    [~, ~, atxbin]= get_perm_prob(x,par);
    [~, ~, atybin]= get_perm_prob(y,par);
    [~, ~, atdxbin]= get_perm_prob(dx,par);
    [~, nbin, atdybin]= get_perm_prob(dy,par);
end
%%
px=zeros(nbin,1); py=px;
for i = 1:nbin
    px(i) = sum(atxbin==i);
    py(i) = sum(atybin==i);
end
px=px/N;
py=py/N;


[Hxy] = jointprob([],[],[],N,atxbin,nbin,atybin,nbin);
[Hdxx] = jointprob([],[],[],N,atdxbin,nbin,atxbin,nbin);
[Hdyy] = jointprob([],[],[],N,atdybin,nbin,atybin,nbin);
Hydxx= jointprob([],[],[],N,atybin,nbin,atdxbin,nbin,atxbin,nbin);
Hxdyy= jointprob([],[],[],N,atxbin,nbin,atdybin,nbin,atybin,nbin);
Hx = getentropy(px);
Hy = getentropy(py);
% Hxy = getentropy(pxy);
% Hdxx = getentropy(pdxx);
% Hdyy = getentropy(pdyy);
% Hydxx = getentropy(pydxx);
% Hxdyy = getentropy(pxdyy);

I12 = Hxy+Hdyy-Hxdyy-Hy;
I21 = Hxy+Hdxx-Hydxx-Hx;
D=(I12-I21)/(I12+I21);
MI = (Hx+Hy-Hxy)/(0.5*(Hx+Hy)); %normalized mutual information to be [0 1] interval

% figure(9)
% 
% pxy=full(pxy);   pxy=reshape(pxy,nbin,nbin);  subplot(2,2,1); surf(pxy); colorbar;title('pxy');
% pxy=full(pdyy);   pxy=reshape(pxy,nbin,nbin);  subplot(2,2,2); surf(pxy); colorbar;title('pdyy');
% pxy=full(pdxx);   pxy=reshape(pxy,nbin,nbin);  subplot(2,2,4); surf(pxy); colorbar;title('pdxx');

