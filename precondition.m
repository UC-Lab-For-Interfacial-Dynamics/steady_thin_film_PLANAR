% boundary condition checker
% 2D planar coordinate system
% Ayaaz Yasin - Nov 13, 2024
% clear; clc; close all;

function Ti = precondition(h,h1,h2)
global F C

%%% set conditions
Tw = C.Tw_func(0);
Ti = 140;
Tv = C.Tv;

alpha = C.a;

fig = figure; fig.Position = [127         458        1113         420];
%%% iterating on Ti
for i = 1:100    

    k   = thermCond(Ti);
    hfg = latentHeat(Ti);
    Pv  = vapPressure(Ti);
    Pd  = disjoinPressure(h);
    [sigma, ~]      = surfTension(Ti);
    kappa           = curvature(h1,h2);
    Pc  = sigma*kappa;
    vl  = specificVol(Ti);

    wayner1 = (Pv*F.M*hfg/(C.Rc*Tv*Ti))*(Ti-Tv);
    wayner2 = (vl*Pv/(C.Rc*Ti))*(Pd + Pc);
    mflux    = (2*C.a/(2-C.a))*sqrt(F.M/(2*pi*C.Rc*Ti))*(wayner1  -  wayner2);
    Ti_star = -(hfg/k)*h*mflux + Tw;    % interface temperature
    
    Ti_new = alpha*Ti_star + (1-alpha)*Ti;
    
    subplot(1,2,1); 
            plot(i, mflux, 'ro', 'LineWidth', 2); hold on; 
            plot(i, wayner1, 'bo', 'LineWidth', 2); hold on; 
            plot(i, wayner2, 'ko', 'LineWidth', 2); hold on; 
            xlabel('iteration'); ylabel('mass flux [kg/m^2-s]');
            set(gca, 'FontWeight','bold','FontSize',14)
    subplot(1,2,2); 
            plot(i, Ti, 'bo', 'LineWidth', 2); hold on; 
            xlabel('iteration'); ylabel('Ti [K]');
            set(gca, 'FontWeight','bold','FontSize',14)
    axis padded; drawnow()
    
    if abs(Ti-Ti_star)<1e-4;     break; end
    Ti = Ti_new;
end


end