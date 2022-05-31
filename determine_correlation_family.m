function [metamodel_pars, aux] = determine_correlation_family(metamodel_pars)

metamodel_pars.correlation_families = { 'Exponential'; 'Gaussian'; 'Linear'; 'Spherical'; 'Cubic'; 'Matern_3_2'; 'Matern_5_2'};


for idx_family = 1:length(metamodel_pars.correlation_families)

    aux(idx_family) = optimize_theta(metamodel_pars, idx_family);
    MLE(idx_family) = aux(idx_family).fit.phi;
%     aux(idx_family).family =  metamodel_pars.correlation_families(idx_family);
end


[~, best] = min(MLE);
metamodel_pars = aux(best);

end

