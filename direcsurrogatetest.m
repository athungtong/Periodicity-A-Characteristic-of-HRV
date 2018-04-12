function [noD uniD12 uniD21 biD D MI hMI]=direcsurrogatetest(I12,I21,I12s,I21s,MI,MIs)

% Test If I12 > I12s at sig lavel =0.05. require number of surrogate = 19.
% Refference: Testing for nonlinearity in time series: the method of surrogate data
% Theiler 1992

noD = 0;   %there is no direction from neither 1->2 nor 2->1
uniD12=0;  %unidirection from 1 to 2
uniD21=0;  %unidirection from 2 to 1
biD=0;   %bidirection but 1->2 is greater than 2->1
D=0; %Directional score of bidirectional 

[h12]=surrcmp(I12,I12s,'g');
[h21]=surrcmp(I21,I21s,'g');

if h12==1 && h21==1
    biD=1;
    D=(I12-I21)/(I12+I21);
elseif h12==1 && h21==0
    uniD12=1; %D=1;
elseif h12==0 && h21==1
    uniD21=1; %D=-1;
elseif h12==0 && h21==0
    noD=1; 
end

if nargin==6
    sigmaMI = (MI-mean(MIs))/std(MIs);
    pMI = erfc(sigmaMI/sqrt(2));
%     hMI=pMI<0.05; %1 if I12>I12s % 0.05 might be a little too much so
%     let's try 0.01
    hMI=pMI<0.01; %1 if MI>MIs 
    if hMI == 0
        MI=0;
    end
else
    MI=NaN;
end