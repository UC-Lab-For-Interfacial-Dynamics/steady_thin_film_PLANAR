% Runge-Kutta Order-4 solver function (multi-dimensional)
% MATH 8010 - Advanced Numerical Analysis
% Ayaaz Yasin - Feb 10, 2023
%
%       [t,y] = ode4(dydt,tSpan,y0,h)
%
% Inputs:
%       dydt ........ function handle
%       tSpan ....... time-span vector 2x1
%       y0 .......... initial condition
%       h ........... time step
% Outputs:
%       t ........... time domain
%       y ........... solution

function [t,y] = ode4(dydt,tSpan,y0,h)

% t = tSpan(1):h:tSpan(2);    % time domain
t = tSpan(1);
y = y0;                     % intialization

for n = 1:length(tSpan(1):h:tSpan(2))-1
    k1 = dydt(t(:,n),            y(:,n),            false);
    k2 = dydt((t(:,n)+h/2),     (y(:,n)+(h/2)*k1),  false);
    k3 = dydt((t(:,n)+h/2),     (y(:,n)+(h/2)*k2),  false);
    k4 = dydt((t(:,n)+h)  ,     (y(:,n)+h*k3),      true);

    y(:,n+1) = y(:,n) + (h/6)*(k1+2*k2+2*k3+k4);
    
    t(n+1) = t(n)+h;

    % figure(1); plot(t(n),y(1,n),'r.'); hold on; drawnow();
    % xlabel('wall distance [m]'); ylabel('film height [m]')

    if isnan(y(1,end));     fprintf('solution is nan!\n');      break; end  % invalid solution
    if y(2,end)>=0;         fprintf('film is flat!\n');         break; end  % flat film
    if ~isreal(y(:,end));   fprintf('solution is not real!\n'); break; end  % solution is complex
    % if y(3,end)<0;          fprintf('concave down!\n');         break; end  % film concave down 

end
end