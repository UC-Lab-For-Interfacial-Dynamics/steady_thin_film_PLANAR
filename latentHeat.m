% Latent Heat of Vaporization Function
% Ayaaz Yasin - Oct 14, 2024
%       hfg = latentHeat(T)

function hfg = latentHeat(T)
    global F
    hfg = polyval(F.hfg_fit, T);  % calculate hfg at given temperature(s)
end