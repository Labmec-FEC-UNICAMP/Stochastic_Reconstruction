function [metamodel_pars] = optimize_theta(metamodel_pars, idx_family)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

metamodel_pars  = start(metamodel_pars, idx_family);
k_max = max(7, min(4, length(metamodel_pars.theta0)));

% This for cannot be paralell because it must respect an execution order,
% i.e., second loop cannot be assessed before first loop.

for i = 1:k_max
    
    metamodel_pars.theta_hat = metamodel_pars.theta;
    metamodel_pars = explore(metamodel_pars, idx_family);
    metamodel_pars = move(metamodel_pars, idx_family);
    metamodel_pars.Delta = circshift(metamodel_pars.Delta,(length(metamodel_pars.Delta)-1));
    
end

end

