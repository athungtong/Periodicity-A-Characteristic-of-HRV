function [Tmax Tmin] = getlocalmaxmin(x)
% [Tmax Tmin] = getlocalmaxmin(x) returns local maxmimum time (Tmax) and
% local minimum time (Tmin) of time series x.

% tmax = zeros(size(x));
% for j = 2:length(x)-1
%     tmax(j) = sign(x(j)-x(j-1)) + sign(x(j)-x(j+1));
% end
% 
% tmax(2:end-1) = sign(x(2:end-1)-x(1:end-2)) + sign(x(2:end-1)-x(3:end));

% % a=sign(x(2:end-1)-x(1:end-2)) ;
% % b=sign(x(2:end-1)-x(3:end));

% temp = sign(diff(x));
% 
% adat=sign(temp(1:end-1));
% bdat = -sign(temp(2:end));

% tmax = -diff(sign(diff(x)));
% Tmax = sign(tmax-2)+1;
% if nargout == 2
%     Tmin = 1-sign(tmax+2); 
% end
[r c]=size(x);


% Tmax = 1 + sign(-diff(sign(diff(x)))-2);
tmax = -diff(sign(diff(x)));

if c==1
    Tmax=[x(1)>x(2);1 + sign(tmax-2); x(end)>x(end-1)];
else
    Tmax=[x(1)>x(2) 1 + sign(tmax-2) x(end)>x(end-1)];
end

if nargout == 2
    if c==1
        Tmin=[x(1)<x(2);1-sign(tmax+2); x(end)<x(end-1)];
    else
        Tmin=[x(1)<x(2) 1-sign(tmax+2) x(end)<x(end-1)];
    end
end

end