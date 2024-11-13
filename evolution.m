% Evolution function
% 2D planar coordinate system
% Ayaaz Yasin - Nov 13, 2024

function dh = evolution(X, H, step)
global F C R
%%% loading and parsing data
% load report.mat
h  = H(1);
h1 = H(2);
h2 = H(3);
h3 = R.h3(end);
Ti_last = R.Ti(end);    % Ti from last full RK step
Tw = C.Tw_func(X);

Ti = Ti_last;
relax = 0.2;
for iter = 1:100     % inner iteration loop
    %%% calling fluid property functions
    rho                 = density(Ti);
    hfg                 = latentHeat(Ti);
    k                   = thermCond(Ti);
    Pv                  = vapPressure(Ti);
    Pd                  = disjoinPressure(h);
    % Pc                  = capPressure(Ti_old, h, h1, h2, h3);
    mu                  = viscosity(Ti);
    [sigma, gamma]      = surfTension(Ti);
    kappa               = curvature(h1,h2);
    Pc                  = sigma*kappa;
    vl                  = specificVol(Ti);
    
    %%% thin film solution
    mflux    = (2*C.a/(2-C.a))*sqrt(F.M/(2*pi*C.Rc*Ti))*((Pv*F.M*hfg/(C.Rc*C.Tv*Ti))*(Ti-C.Tv)  -  (vl*Pv/(C.Rc*Ti))*(Pd + Pc));      % mass flux. CHANGE TO BETA-FORM
    Ti_star       = -(hfg/k)*h*mflux + Tw;    % interface temperature
    Ti_new = relax*Ti_star + (1-relax)*Ti;

    if abs(Ti-Ti_star)<1e-4;     break; end
    Ti = Ti_new;
end
    % if Ti < 0; error('interface temp, Ti, is negative!'); end
dT_dx    = (Ti-Ti_last)/C.dx;   % temperature gradient
dPd_dx   = -3*F.A*h1/h^4;     % disjoining pressure gradient
ds_dx    = gamma*dT_dx;     % surface tension gradient


dPl_dx = 6*(mflux/h)*(mu/rho) - 3*h*ds_dx;   % liquid pressure gradient

% film height derviatives
dh = nan(3,1); 
dh(1) = H(2);
dh(2) = H(3);
dh(3) = 3*(dh(2)^2)*dh(1)/(1+(dh(1))^2) -(gamma/sigma)*dT_dx*dh(2) - (1+(dh(1))^2)^(3/2)*(1/sigma)*(dPl_dx + dPd_dx);

if step % step is passed by the ode solver. it is true if this call to the evolution function is NOT at an intermediate RK-step.
    %%% save external data
    R.n(end+1)         = R.n(end)+1;
    R.xEval(end+1)     = X;
    R.h3(end+1)        = dh(3);
    R.Ti(end+1)        = Ti;
    R.Tw(end+1)        = Tw;
    R.mflux(end+1)     = mflux;
    R.kappa(end+1)     = kappa;
    R.Pc(end+1)        = Pc;
    R.Pd(end+1)        = Pd;
    
    %%% progress update
    fprintf('%i\t\t%0.2f\t\t%0.16e\n', R.n(end-1), (X-C.x_in)/C.L, h)
end

end