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
%       varargin{1} ..... optional input: initial guess for interface temp
%                                  if no optional input, (Tw+Tv)/2 is used.
% Outputs:
%       Ti .............. interface temperature [K]
%       mflux ........... phase change mass flux [kg/m^3]


function [Ti, mflux] = precondition(h,h1,h2, Tv, Tw, varargin)
global F C          % load global data variables, F-fluid data, C-constants
if nargin==5
    Ti = (Tv+Tw)/2;     % initial guess for interface temperature Ti
elseif nargin==6
    Ti = varargin{1};
else
    error('Precondition function: Too many inputs.')
end
relax = 1e-1;        % relaxation factor for temperature
maxIter = 1e4;      % max iterations allowed
minIter = 30;

for i = 1:maxIter    
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
    waynerCoeff = (2*a/(2-a))*sqrt(F.M/(2*pi*C.Rc*Ti));
    term1  = waynerCoeff*(Pv_sat*F.M*hfg/(C.Rc*Tv*Ti))*(Ti-Tv); % wayner thermal term
    term2  = waynerCoeff*(vl*Pv_sat/(C.Rc*Ti))*(Pd + Pc);       % wayner pressure term
    mflux  = term1  -  term2;  % mass flux

    %%% temperature
    Ti_star = -(hfg/k)*h*mflux + Tw;             % interface temperature
    Ti_new  = relax*Ti_star + (1-relax)*Ti;       % relaxation

    if strcmpi(C.debug,'on') % plots only if debug mode is 'on'
        fig = figure(5); fig.Position = [127,458,1113,420];
        subplot(1,2,1); 
                plot(i, mflux, 'ko', 'LineWidth', 2); hold on; 
                plot(i, term1, 'ro', 'LineWidth', 2); hold on; 
                plot(i, term2, 'bo', 'LineWidth', 2); hold on; 
                xlabel('precondition iteration'); ylabel('mass flux [kg/m^2-s]');
                set(gca, 'FontWeight','bold','FontSize',14)
                % legend('wayner','term1','term2','location','best')
        subplot(1,2,2); 
                plot(i, Ti, 'ro', 'LineWidth', 2); hold on; 
                % yline(Tv,'b--','LineWidth',2);
                % yline(Tw,'k--','LineWidth',2);
                xlabel('precondition iteration'); ylabel('Ti [K]');
                set(gca, 'FontWeight','bold','FontSize',14)
                % legend('interface','vapor','wall','location','best')
        axis padded; hold on;drawnow();
    end

    if abs(Ti-Ti_star)<1e-4 && i>minIter;     break; end     % stopping criterion
    Ti = Ti_new;                                % variable shift for next iteration
    if i==maxIter; error('warning: precondition loop maxed out. increase max iterations.'); end
end
if strcmpi(C.debug,'on')
    input('Press ''Enter'' to continue...','s'); 
    close(5);
end
end