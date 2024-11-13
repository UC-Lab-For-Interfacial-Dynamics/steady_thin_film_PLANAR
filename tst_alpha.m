% Transition state theory based accommodation coefficient calculator
% Ayaaz Yasin - Oct 31, 2022
% Inputs:
% rv...........vapor density
% rl...........liquid density
% Outputs:
% a............accommodation coefficient
% l............non-dimensional length-scale
% Reference: nagayama_2003
function [a,l] = tst_alpha(rg, rl)
    l    = nthroot(rg./rl,3);
    linv = nthroot(rl./rg,3);
    a    = (1-l).*exp((-1/2).*1./(linv-1));
end