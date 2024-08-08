function out = defaultVaporsTransportEquation(vaporConc,dVaporConc,param)

% 'out' is the output of the function that denotes the d vector of the
% vapor concentrations after the function. If 'out' equals 'dVaporConc', 
% this function has no effect because 'dVaporConc' is the d vector of the
% vapor concentrations before the function.

% initialize output vector
out = dVaporConc;

% example where vapor 1 turns into vapor 2 in one hour
out(1) = out(1) - vaporConc(1)/3600;
out(2) = out(2) + vaporConc(1)/3600;

