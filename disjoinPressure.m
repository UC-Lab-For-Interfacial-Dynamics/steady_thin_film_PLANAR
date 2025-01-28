% Disjoining Pressure Function
%       Calculates the local disjoining pressure as the function of the
%       film height, h. Hamaker constant is needed from the global 
%       variable 'F.A'. The exponent is needed from the global 
%       variable 'C.pdexp'.
%       
% Ayaaz Yasin - Sep 11, 2024
% 
%       Pd = disjoinPressure(h)

function Pd = disjoinPressure(h)
    global F C
    Pd = F.A/(h^3);           % disjoining pressure
end