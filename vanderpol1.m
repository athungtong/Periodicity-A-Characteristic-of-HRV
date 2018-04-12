function [y time]=vanderpol1(Mu,w,d)
%%

    tau=0.2;     Maxtime=100;
    time = (0:tau:Maxtime-tau)';
    noise = mvnrnd(0,1,length(time));
     %noise = rand(length(time));

    y0 = [2; 0];

    ode = @(t,y) vanderpoldemo(t,y);
 
    [t,y] = ode45(ode, time, y0);

    % Plot of the solution
%    figure()
%    plot(t,y(:,1))
%     xlabel('t')
%     ylabel('solution y')
%     title('van der Pol Equation, \mu = 1')
    function dydt = vanderpoldemo(t,y)
        %VANDERPOLDEMO Defines the van der Pol equation for ODEDEMO.

        % Copyright 1984-2002 The MathWorks, Inc.
        % $Revision: 1.2 $  $Date: 2002/06/17 13:20:38 $
         dydt = [y(2); Mu*(1-y(1)^2)*y(2)-(w+d*noise(round(t/tau)+1))*y(1)];
    end
end

    
