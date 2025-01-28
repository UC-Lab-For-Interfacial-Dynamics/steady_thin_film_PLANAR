% Precondition function
%       Function to calculate the interface temperature and mass flux at
%       any point along the film. The mass flux and 1D heat conduction is
%       solved iteratively. Due to the non-linearity, the temperature
%       calculation is under-relaxed.
% Ayaaz Yasin - Nov 13, 2024
%
%       [Ti, mflux] = precondition(h,h1,h2, Tv, Tw)
%
% Inputs:
%       h, h1, h2 ....... film height & derivatives, h, h', h'' [m, -, 1/m]
%       Tv .............. vapor temperature [K]
%       Tw .............. wall temperature [K]
% Outputs:
%       Ti .............. interface temperature [K]
%       mflux ........... phase change mass flux [kg/m^3]


function [Ti, mflux] = precondition(h,h1,h2, Tv, Tw)
global F C          % load global data variables, F-fluid data, C-constants
Ti = (Tv+Tw)/2;     % initial guess for interface temperature Ti
relax = 0.1;        % relaxation factor for temperature

for i = 1:100    
    %%% fluid properties at the current Ti
    k           = thermCond(Ti);
    hfg         = latentHeat(Ti);
    Pv_sat      = vapPressure(Ti);  % saturation pressure from NIST curvefit
    Pd          = disjoinPressure(h);
    [sigma, ~]  = surfTension(Ti);
    kappa       = curvature(h1,h2);
    Pc          = sigma*kappa;
    vl          = specificVol(Ti);
    [rhov,rhol] = density(Ti);
    [a,~]       = tst_alpha(rhov,rhol);
    
    %%% mass flux
    term1  = (Pv_sat*F.M*hfg/(C.Rc*Tv*Ti))*(Ti-Tv); % wayner thermal term
    term2  = (vl*Pv_sat/(C.Rc*Ti))*(Pd + Pc);       % wayner pressure term
    mflux  = (2*a/(2-a))*sqrt(F.M/(2*pi*C.Rc*Ti))*(term1  -  term2);  % mass flux

    %%% temperature
    Ti_star = -(hfg/k)*h*mflux + Tw;             % interface temperature
    Ti_new  = relax*Ti_star + (1-relax)*Ti;       % relaxation

    if strcmpi(C.debug,'on') % plots only if debug mode is 'on'
        fig = figure(5); fig.Position = [127,458,1113,420];
        subplot(1,2,1); 
                plot(i, mflux, 'ro', 'LineWidth', 2); hold on; 
                xlabel('precondition iteration'); ylabel('mass flux [kg/m^2-s]');
                set(gca, 'FontWeight','bold','FontSize',14)
        subplot(1,2,2); 
                plot(i, Ti, 'ro', 'LineWidth', 2); hold on; 
                xlabel('precondition iteration'); ylabel('Ti [K]');
                set(gca, 'FontWeight','bold','FontSize',14)
        axis padded; drawnow()
    end

    if abs(Ti-Ti_star)<1e-4;     break; end     % stopping criterion
    Ti = Ti_new;                                % variable shift for next iteration
end
try; close(5); end
end