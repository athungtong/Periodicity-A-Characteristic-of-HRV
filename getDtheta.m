function [theta1,theta2,Dtheta1,Dtheta2]=getDtheta(theta1,theta2,nmseries,safr,tau)
if nargin<5
    tau=1/safr;
end
% if size(nmseries,2)>1 
% %    theta1 = (nmseries(:,1).*theta1);    theta2 = (nmseries(:,2).*theta2);
% %    theta1 = theta1./nmseries(:,2).*nmseries(:,1);    theta2 = theta2./nmseries(:,1).*nmseries(:,2);
% end
% plot(theta1,'.-');hold on; plot(theta2,'r.-');hold off
% stop
N=round(tau*safr); if N==0, N=1;end
Dtheta1=(theta1(1+N:end)-theta1(1:end-N))/tau;%.*nmseries(1+N:end,1);
Dtheta2=(theta2(1+N:end)-theta2(1:end-N))/tau;%.*nmseries(1+N:end,2);
theta1=theta1(1:end-N);  theta2=theta2(1:end-N);

% use theta(t+tau) instead of incremence of theta
% Dtheta1=theta1(1+N:end);%.*nmseries(1+N:end,1);
% Dtheta2=theta2(1+N:end);%.*nmseries(1+N:end,2);
% theta1=theta1(1:end-N);  theta2=theta2(1:end-N);


