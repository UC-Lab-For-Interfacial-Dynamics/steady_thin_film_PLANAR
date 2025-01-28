% Runge-Kutta Order-6 solver function (multi-dimensional)
%       Fixed-step 6th-order Runge-Kutta solver. This function can handle
%       ODEs of order>1. If dydt=y', then y=y. 
%       If dydt=[y',y'',y'''], then y=[y,y',y''].
%       Unlike a standard ODE solver, this function:
%           (1) Calls the derivative function with three inputs
%               dydt(t, y, step) where step is 'false' for k1, k2, k3, but
%               'true' for k4. This lets the derivative function know when
%               it is being called at actual solution step and data can be
%               saved for post-processing.
%           (2) Stopping conditions check for invalid and complex solutions
%               as well as, the first derivative becoming zero or positive.
% Ayaaz Yasin - Oct 30, 2024
%
%       [t,y] = ode6(dydt,tSpan,y0,h)
%
% Inputs:
%       dydt ........ function handle
%       tSpan ....... time-span vector 2x1
%       y0 .......... initial condition
%       h ........... time step
% Outputs:
%       t ........... time domain
%       y ........... solution

function [t, y] = ode6(dydt, tSpan, y0, h)
    t = tSpan(1);
    y = y0;

    for n = 1:length(tSpan(1):h:tSpan(2))-1
        % Intermediate steps for the RK solver
        k1 = dydt(t(:,n),             y(:,n),                       false);
        k2 = dydt(t(:,n) + h/3,       y(:,n) + h/3*k1,              false);
        k3 = dydt(t(:,n) + h/3,       y(:,n) + h/6*k1   + h/6*k2,   false);
        k4 = dydt(t(:,n) + h/2,       y(:,n) + h/8*k1   + 3*h/8*k3, false);
        k5 = dydt(t(:,n) + 2*h/3,     y(:,n) - 3*h/8*k2 + 9*h/8*k4, false);
        k6 = dydt(t(:,n) + h,         y(:,n) + h*(k1 - 4*k3 + 4*k4 - k5), true);

        % Integration
        y(:,n+1) = y(:,n) + h*(k1 + 4*k3 + 4*k4 + k5 + k6) / 15;
        t(n+1) = t(n) + h;

        % Stopping Conditions
        if isnan(y(1,end))                      % invalid solution
            fprintf('Solution terminated!\nFilm height became invalid (nan) at x = %0.2f nm', t(end)*1e9);      
            break; 
        end  
        if ~isreal(y(:,end))                    % solution is complex
            fprintf('Solution terminated!\nFilm height or derivatives became non-real at x = %0.2f nm\n', t(end)*1e9); 
            break; 
        end  
        if y(2,end)>=0                          % film is flat
            fprintf('Solution terminated! \nFilm became flat at x = %0.2f nm.\n', t(end)*1e9);         
            break; 
        end  
    end
end
