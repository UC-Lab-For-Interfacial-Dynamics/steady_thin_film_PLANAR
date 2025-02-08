% Disjoining Pressure Function
%       Calculates the local disjoining pressure as the function of the
%       film height, h. Hamaker constant is needed from the global 
%       variable 'F.A'.
%       
% Ayaaz Yasin - Sep 11, 2024
% 
%       Pd = disjoinPressure(h)

function Pd = disjoinPressure(h)
    global F
    Pd = F.A./(h.^3);           % disjoining pressure
end