function yhat=mypolyfit(y,d,T,A,N)
% yhat=mypolyfit(y,T,A,D,N) return estimated yhat which is the polynomial fitted 
% version of input y at time T (vector). A is the matrix, D is the degree
% of polynomial and N is the length of the data

switch nargin
    case{1}
        d=1; [T A N] = genAmatrix(y,d);
    case{2}
        [T A N] = genAmatrix(y,d);
    case{4}
        N=length(y);
end

yhat=zeros(N,d);
d0=d;
if length(y)<=d
   d = length(y)-1; [T A N] = genAmatrix(y,d);
end

if d==0,yhat(:,1:d)=y;return; end

phi = zeros(d+1,1);
b=[sum(y)/N;(y'*T)'/N];

%% For D <4, we can solve equation analytically
% For d=1
phi(1)=b(1);
phi(2)=b(2)/A(1);
yhat(:,1)=phi(1)+T(:,1)*phi(2);
if d>1% For d=2
    %phi(2)=previous
    phi(3) = (b(3)-A(1)*b(1))/(A(2)-A(1)^2);
    phi(1) = b(1)-A(1)*phi(3);
    yhat(:,2)=phi(1)+T(:,1:2)*phi(2:3);
    if d>2
        %For d=3;
        % phi(1) and phi(3) are equal to previous
        phi(4) = (A(2)*b(2)-A(1)*b(4))/(A(2)^2-A(1)*A(3));
        phi(2) = (b(2)-A(2)*phi(4))/A(1);
        yhat(:,3)=phi(1)+T(:,1:3)*phi(2:4);
        if d>3
            %% This code is for D>3
            %compute b matrix for many degree
            %the code is not done, we need to compute A in the different way. 
            %Need to correct the code in getZmatrix
            yhat=zeros(N,d);
            b=zeros(d+1,1);

            b(1)=sum(y)/N;
            for d=1:d
                b(d+1)=sum(y.*T(:,d))/N; % sum(y.*t.^d);
                row=1:d+1;
                [Q R] = qr(A(row,row));
                bb = Q'*b(row);    
                phi = zeros(length(R),1); %This is the kernal
                phi(end) = bb(end)/R(end,end);
                for i=length(R)-1:-1:1
                   phi(i) =  (bb(i) - R(i,i+1:end)*phi(i+1:end) )/R(i,i);
                end  
                yhat(:,d)=phi(1)+T(:,1:d)*phi(2:end);    %polyval(flipud(phi),t);
            end
        end
    end
end
if d0>d
    yhat(:,d+1:d0)=yhat(:,d);
end

end

function [T A N] = genAmatrix(y,D)
    N=length(y);
    t=cell(2*D,1);
    t{1} = (-1+1/N:2/N:1-1/N)';
    for i=2:2*D
        t{i}=t{i-1}.*t{1};
    end
    
    T=[]; A=[];
    for i=1:D
        A = [A;sum(t{i*2})/N];
        T=[T t{i}];
    end
end
