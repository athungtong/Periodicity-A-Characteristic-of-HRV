function [PoinCare]=PoincarePlotCalculation(X,tau)
% compute poincare of X with time delay tau
[SD1C SD2C SD_Ratio CDOWN CUP]=PoinCareFnct(X,tau);
PoinCare.SD1=SD1C;
PoinCare.SD2=SD2C;
PoinCare.SDRatio=SD_Ratio;
PoinCare.UP=CUP;
PoinCare.DOWN=CDOWN;

end

function [SD1C SD2C SD_Ratio CDOWN CUP]=PoinCareFnct(X,tau)
%%
x = X(1:end-tau);  y = X(1+tau:end);

% Step 1. Find the new coordinate which is -45 degree from the original coordinate
xcorr = [1; 1]/sqrt(2); ycorr = [-1; 1]/sqrt(2);
% Step 2. Change coordinate of x and y
xdat = [x y]*xcorr; ydat = [x y]*ycorr;
% Step 2.5 Remove artifact
% Artifact is defined as any point ofwhich Cook's distance>4/n
% Where n is the number of point 
% Reference Bollen, Kenneth A.; and Jackman, Robert W. (1990); 
% Regression diagnostics: An expository treatment of outliers 
% and influential cases, in Fox, John; and Long, J. Scott (eds.); 
% Modern Methods of Data Analysis (pp. 257-91). Newbury Park, CA: Sage
% -- Editted 6/7/2014 By Anurak Thungtong
hii = leverage(xdat);   ri = xdat-mean(xdat);   ci = hii./(1-hii).*ri.^2;
indx=ci>4/length(xdat); % index of outliers

hii = leverage(ydat);   ri = ydat-mean(ydat);   ci = hii./(1-hii).*ri.^2;
indy=ci>4/length(ydat); % index of outliers

index=unique([find(indx);find(indy)]);
xdat(index)=[]; 
ydat(index)=[]; 

% Step 3 and 4. SD1 is the SD of ydat and SD2 is the SD of xdat
SD1C = std(ydat);   SD2C = std(xdat);
% Step 5. Compute SD_Ratio
if SD2C==0
    SD_Ratio=Inf;
else
    SD_Ratio=SD1C/SD2C;
end
% Step 6. Check symmetry around xcorr axis
ydat_up = ydat(ydat>=0);
ydat_down = ydat(ydat<0);
CUP = sum(ydat_up.^2)/sum(ydat.^2);
CDOWN = sum(ydat_down.^2)/sum(ydat.^2);
% 
% figure(3)
% % % subplot(133);
% plot(x,y,'.')
% hold on
% line([min(x) max(x)],[min(y) max(y)],'color','r');hold off
% xlabel('x(t)');
% ylabel('x(t+1)');
% % con=input('con?');

% % title(label);
% out=[SD1C SD2C SD_Ratio];
%     out
% 
% s=input('save?');
% if s
%     xy=[x y];
%     xlswrite('poincare.xls',xy,'bsl');
%     'optimal SD'
%     out
% end
%This is the same as above but more difficult to understand
% L = length(x);
% SD1C = sqrt((1/L)*sum(((x-y)-mean(x-y)).^2)/2);
% SD2C = sqrt((1/L)*sum(((x + y) - mean(x + y)).^2)/2);
% SD1I = sqrt((1/L) * (sum((x - y).^2)/2));
% xy = (y-x)/sqrt(2);
% indices_up = xy > 0;
% indices_down = xy < 0;
% SD1UP = sqrt(sum(xy(indices_up).^2)/L);
% SD1DOWN = sqrt(sum(xy(indices_down).^2)/L);
% CUP = SD1UP^2/SD1I^2;
% CDOWN = SD1DOWN^2/SD1I^2;
% SD_Ratio=SD1C/SD2C;
end
