function [ metamodel_pars ] = explore(metamodel_pars, idx_family)

V = metamodel_pars.theta;

for i = 1:length(metamodel_pars.theta)
    
    if metamodel_pars.theta(i) == metamodel_pars.theta_lower(i)
        V(i) = metamodel_pars.theta_lower(i)*sqrt(metamodel_pars.Delta(i));
        flag = true;
        
    elseif metamodel_pars.theta(i) == metamodel_pars.theta_upper(i)
        V(i) = metamodel_pars.theta_upper(i)/sqrt((metamodel_pars.Delta(i)));
        flag = true;
        
    else
        V(i) = min(metamodel_pars.theta(i)*metamodel_pars.Delta(i), metamodel_pars.theta_upper(i));
        flag = false;
    end
    
    metamodel_pars_V = metamodel_pars;
    metamodel_pars_V.theta = V;
    metamodel_pars = objective_function( metamodel_pars, idx_family );
    metamodel_pars_V = objective_function( metamodel_pars_V, idx_family );
    
    if metamodel_pars_V.fit.phi < metamodel_pars.fit.phi
        metamodel_pars.theta = metamodel_pars_V.theta;
    else
        if flag == false
            
            V(i) = max(metamodel_pars.theta(i)/metamodel_pars.Delta(i), metamodel_pars.theta_lower(i));
            metamodel_pars_V.theta = V;
            metamodel_pars = objective_function( metamodel_pars, idx_family );
            metamodel_pars_V = objective_function( metamodel_pars_V, idx_family );
            
            if metamodel_pars_V.fit.phi < metamodel_pars.fit.phi
                metamodel_pars.theta = metamodel_pars_V.theta;
            end
        end
    end
end

end

