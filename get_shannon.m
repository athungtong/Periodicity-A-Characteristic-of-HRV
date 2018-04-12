function [IXY HX HY HXY] = get_shannon(px,py,pxy,normal)
    
    if nargin == 3, normal = 1;end
    
    [nxbin nybin] = size(pxy);    
    px(px==0)=[];
    py(py==0)=[];
    pxy=reshape(pxy,nxbin*nybin,1);
    pxy(pxy==0)=[];
    
    if normal == 0 || normal == 1
        HX=-dot(px,log2(px));    
        HY=-dot(py,log2(py));
        HXY=-dot(pxy,log2(pxy));
        if normal == 0
            % no normalization required. IXY = 0 up to some positive number
            IXY = (HX+HY-HXY);
        else
            % normalized mutual information to be [0 1] interval
            % ref: Zhi-Yuan Su; Tzuyin Wu, Yeng-Tseng Wang, Hsin-Yi Huang. An
            % investigation into the linear and non linear correlation of two music
            % walk sequences
            % Physica D 237 (2008) 1815-1824
            % IXY_norm = (HX+HY-HXY)/(0.5*(HX+HY));
            IXY = (HX+HY-HXY)/(0.5*(HX+HY));
        end
    else
        % normalized mutual information to be [0 1] interval method 2
        % basic idea is to use log base equal to number of bin
        % Q. Wang et al. Physica D 200 (2005) 287-295
        %important note: not gaurantee mi=1 for the identical signal
        HX=-dot(px,(log(px)/log(nxbin))); %Change log base to the # of bin
        HY=-dot(py,(log(py)/log(nybin)));
        HXY=-dot(pxy,(log(pxy)/log(nxbin)));    
        IXY = (HX+HY-HXY);
    end
end

% % normalized mutual information to be [0 1] interval
% % ref: Zhi-Yuan Su; Tzuyin Wu, Yeng-Tseng Wang, Hsin-Yi Huang. An
% % investigation into the linear and non linear correlation of two music
% % walk sequences
% % Physica D 237 (2008) 1815-1824
% % IXY_norm = (HX+HY-HXY)/(0.5*(HX+HY));
% function [IXY HX HY HXY] = get_shannon_norm1(px,py,pxy)
%     [nxbin nybin] = size(pxy);    
%     px(px==0)=[];
%     HX=-dot(px,log2(px));
%     py(py==0)=[];
%     HY=-dot(py,log2(py));
%     pxy=reshape(pxy,nxbin*nybin,1);
%     pxy(pxy==0)=[];
%     HXY=-dot(pxy,log2(pxy));    
%     IXY = (HX+HY-HXY)/(0.5*(HX+HY));
%         
% end
% 
% 
% % normalized mutual information to be [0 1] interval method 2
% % basic idea is to use log base equal to number of bin
% % Q. Wang et al. Physica D 200 (2005) 287-295
% %important note: not gaurantee mi=1 for the identical signal
% function [IXY HX HY HXY] = get_shannon_norm2(px,py,pxy)
%     [nxbin nybin] = size(pxy);    
%     px(px==0)=[];
%     HX=-dot(px,(log(px)/log(nxbin))); %Change log base to the # of bin
%     py(py==0)=[];
%     HY=-dot(py,(log(py)/log(nybin)));
%     pxy=reshape(pxy,nxbin*nybin,1);
%     pxy(pxy==0)=[];
%     HXY=-dot(pxy,(log(pxy)/log(nxbin)));
%     
%     IXY = (HX+HY-HXY);
%         
% end