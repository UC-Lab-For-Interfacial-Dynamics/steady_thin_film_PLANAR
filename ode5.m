% Runge-Kutta Order-5 solver function (multi-dimensional)
% Fixed-Step solver
% Ayaaz Yasin - Oct 30, 2024
%
%       [t,y] = ode5(dydt,tSpan,y0,h)
%
% Inputs:
%       dydt ........ function handle
%       tSpan ....... time-span vector 2x1
%       y0 .......... initial condition
%       h ........... time step
% Outputs:
%       t ........... time domain
%       y ........... solution

function [t, y] = ode5(dydt, tSpan, y0, h)
    t = tSpan(1);
    y = y0;

    for n = 1:length(tSpan(1):h:tSpan(2)) - 1
        k1 = dydt(t(:, n),            y(:, n),            false);
        k2 = dydt(t(:, n) + h/4,      y(:, n) + (h/4)*k1, false);
        k3 = dydt(t(:, n) + h/4,      y(:, n) + (h/8)*k1    + (h/8)*k2, false);
        k4 = dydt(t(:, n) + h/2,      y(:, n) - (h/2)*k2    + (h)*k3, false);
        k5 = dydt(t(:, n) + 3*h/4,    y(:, n) + (3*h/16)*k1 + (9*h/16)*k4, false);
        k6 = dydt(t(:, n) + h,        y(:, n) - (3*h/7)*k1  + (2*h/7)*k2 + 12*h/7*k3 - 12*h/7*k4 + 8*h/7*k5, true);

        y(:, n+1) = y(:, n) + (7*k1 + 32*k3 + 12*k4 + 32*k5 + 7*k6)*(h/90);
        t(n+1) = t(n) + h;

        % figure(1); plot(t(n),y(1,n),'r.'); hold on; drawnow();
        % xlabel('wall distance [m]'); ylabel('film height [m]')
    
        if isnan(y(1,end));     fprintf('solution is nan!\n');      break; end  % invalid solution
        if y(2,end)>=0;         fprintf('film is flat!\n');         break; end  % flat film
        if ~isreal(y(:,end));   fprintf('solution is not real!\n'); break; end  % solution is complex
        if y(3,end)<0;          fprintf('concave down!\n');         break; end  % film concave down 
    end
end
