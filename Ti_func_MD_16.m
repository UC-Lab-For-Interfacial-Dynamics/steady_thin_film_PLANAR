% Interface temperature function
%       Polynomial fit of the interface tempeature from MD.
% Ayaaz Yasin - Apr 18, 2025
%
%       Ti = Ti_func_MD(x)
%
% Input:    x_in .................. location along the interface [meters]
% Output:   Ti .................... vapor temperature [K]

function Ti = Ti_func_MD_16(X)
    X = X*1e9 + 14.4;   % offsetting input and converting to nm
    Ti = [];
    for x=X
    if 14.3<=x && x<=16.5
        coeff = [9.76832482e-01, -4.34173931e+01, 6.44884732e+02, -3.07596085e+03];
    elseif 16.5<x && x<=32.0852
        coeff = [4.75816119e-04, -5.97219668e-02, 2.97053708e+00, -7.32473687e+01, 8.96796801e+02, -4.22241621e+03];
    else 
        coeff = nan;
    end
    T = polyval(coeff, x);
    Ti(end+1)=T;
    end

end