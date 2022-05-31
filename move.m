function [ metamodel_pars ] = move(metamodel_pars, idx_family)

metamodel_pars_V = metamodel_pars;
metamodel_pars = objective_function( metamodel_pars, idx_family );
v = metamodel_pars.theta./metamodel_pars.theta_hat;
if all(v==1)
    metamodel_pars.Delta = metamodel_pars.Delta.^(1/5);
end

rpt = 1;
while rpt
    V =  min(max(metamodel_pars.theta.*v, metamodel_pars.theta_lower), metamodel_pars.theta_upper);
    metamodel_pars_V.theta = V;
    metamodel_pars_V = objective_function( metamodel_pars_V, idx_family );
    if metamodel_pars_V.fit.phi < metamodel_pars.fit.phi
        metamodel_pars.theta = metamodel_pars_V.theta;
        v = (v.^2);
        metamodel_pars.fit =  metamodel_pars_V.fit;
    else
        rpt = 0;
    end
    if any((V ==  metamodel_pars.theta_lower) | (V ==  metamodel_pars.theta_upper))
        rpt = 0;
    end 
end

metamodel_pars.Delta = (metamodel_pars.Delta).^(1/4);

end

