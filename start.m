function [ metamodel_pars ] = start(metamodel_pars, idx_family)

q = length(metamodel_pars.theta0);
theta = metamodel_pars.theta0;
Delta = transpose(2.^((1:1:q)./(q+2)));
pos = find((metamodel_pars.theta_lower - metamodel_pars.theta_upper) == 0);
if ~isempty(pos)
Delta(pos) = 1;
theta(pos) = upper(pos);
clear pos
end
pos = find(theta < metamodel_pars.theta_lower | theta > metamodel_pars.theta_upper);
theta(pos) = (metamodel_pars.theta_lower(pos).*metamodel_pars.theta_upper(pos).^7).^(1/8);
metamodel_pars.theta = theta;

metamodel_pars = objective_function( metamodel_pars, idx_family );

if length(pos)>1
    flag = [];
    theta_hat = theta; phi_hat = metamodel_pars.fit.phi;
    J = pos(1);
    
    for i = 1:length(pos)
        j = pos(i);
        phi_bar = phi_hat;
        v = ones(length(pos),1)*0.5;
        v(j) = 1/16;
        alpha = min(log(metamodel_pars.theta_lower(pos)./theta_hat(pos)) ./log(v'));
        v = transpose(v.^(alpha/5));
        for k = 1:4
            V = (v.^k).*theta_hat;
            metamodel_pars_V = metamodel_pars;
            metamodel_pars_V.theta = V;
            metamodel_pars_V = objective_function( metamodel_pars_V, idx_family );
            if metamodel_pars_V.fit.phi < phi_bar
                metamodel_pars_V.fit.phi = phi_bar;
                
                if metamodel_pars_V.fit.phi < metamodel_pars.fit.phi
                    theta = V;
                    J = j;
                else
                    flag = 0;
                end
            end
        end
    end
    if ~isempty(flag)
        Delta([1 J]) = Delta([J 1]);
    end
end

metamodel_pars.theta = theta;
metamodel_pars.Delta = Delta;

end

