% Evolution function
%       This function serves as the dervivative function for the ODE
%       solution. It calculates the film derivatives 'dh' at the current
%       step using data from the previous step. Data is saved to the global
%       report variable 'R' only if the current call is from the final RK
%       step, i.e. step='true'.
% Ayaaz Yasin - Nov 13, 2024
%
%       dh = evolution(X, H, step)
%
% Inputs:
%       X ................. x-location [m]
%       H ................. film height vector, [h, h', h'']
% Output:
%       dh ................ film height derivative vector, [h', h'', h''']

function dh = evolution(X, H, step)

%%% Loading previous step data
    global F C R            % global variables F-fluid data, C-constants, R-reports
    h  = H(1);              % first derivative
    h1 = H(2);              % second derivative
    h2 = H(3);              % third derivative
    Ti_last = R.Ti(end);    % interface temperature from last full RK step (k4)
    Tw = C.Tw_func(X);      % wall temperature at the current step
    Tv = C.Tv_func(X);      % vapor temperature at the current step

%%% Inner loop to solve for interface temperature and mass flux
    [Ti, mflux] = precondition(h,h1,h2,Tv,Tw);

%%% Evolution equation
    % fluid properties
    [rhov,rhol]            = density(Ti);
    Pd                  = disjoinPressure(h);
    mu                  = viscosity(Ti);
    [sigma, gamma]      = surfTension(Ti);
    kappa               = curvature(h1,h2);
    Pc                  = sigma*kappa;
    a                   = tst_alpha(rhov, rhol);

    % derivatives 
    dT_dx    = (Ti-Ti_last)/C.dx;                   % temperature gradient
    dPd_dx   = -3*F.A*h1/h^4;                       % disjoining pressure gradient
    ds_dx    = gamma*dT_dx;                         % surface tension gradient
    dPl_dx   = 6*(mflux/h)*(mu/rhol) - 3*h*ds_dx;   % liquid pressure gradient

    % film height derviatives from the evolution equation
    dh = nan(3,1); 
    dh(1) = H(2);
    dh(2) = H(3);
    dh(3) = 3*(dh(2)^2)*dh(1)/(1+(dh(1))^2) -(gamma/sigma)*dT_dx*dh(2) - (1+(dh(1))^2)^(3/2)*(1/sigma)*(dPl_dx + dPd_dx); % evolution equation

%%% Save solution data
    if step % only saveed if step=true, i.e. this call to the evolution function is NOT at an intermediate RK-step.
        % save data to global report variable 'R'
        R.n(end+1)         = R.n(end)+1;
        R.xEval(end+1)     = X;
        R.h3(end+1)        = dh(3);
        R.Ti(end+1)        = Ti;
        R.Tw(end+1)        = Tw;
        R.Tv(end+1)        = Tv;
        R.mflux(end+1)     = mflux;
        R.kappa(end+1)     = kappa;
        R.Pc(end+1)        = Pc;
        R.Pd(end+1)        = Pd;
        R.a(end+1)         = a;

        % print progress to command window (if R.print=true)
        if strcmpi(C.debug,'on')    % print film height in debug mode
            fprintf('%i\t\t%0.2f\t\t%0.16e\n', R.n(end-1), (X-C.x_in)/C.L, h); 
        else                        % else print status bar
            status = floor((X-C.x_in)/C.L*50);
            bar = strcat('[',repmat('#', 1, status), repmat('-', 1, 50-status),']');
            clc; fprintf('%s\n',bar)
        end
    end
end