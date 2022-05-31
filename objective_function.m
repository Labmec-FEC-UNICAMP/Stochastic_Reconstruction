function [ metamodel_pars ] = objective_function( metamodel_pars, idx_family )

%% Avoiding broadcasts variables aiming to use parfor

scaled_X = metamodel_pars.scaled_X;
scaled_Y = metamodel_pars.scaled_Y;
theta = metamodel_pars.theta;
theta0 = metamodel_pars.theta0;


%%
switch metamodel_pars.correlation_families{idx_family}
    case 'Exponential'
        K = @(theta, distance_matrix) (exp(-distance_matrix./theta));
    case 'Gaussian'
        K = @(theta, distance_matrix) (exp(-distance_matrix.^2./theta.^2));
    case 'Linear'
        K = @(theta, distance_matrix) (max(0, 1 - distance_matrix./theta));
    case 'Spherical'
        %         qsi = min(1, theta*distance_matrix);
        K = @(theta, distance_matrix) 1 - 1.5.*(min(1, distance_matrix./theta)) + 0.5*(min(1, distance_matrix./theta)).^3;
    case 'Cubic'
        %         qsi = min(1, theta*distance_matrix);
        K = @(theta, distance_matrix) 1 - 3.*(min(1, distance_matrix./theta)).^2 + 2*(min(1, distance_matrix./theta)).^3;
    case 'Matern_3_2'
        K = @(theta, distance_matrix)...
            (1 + sqrt(3).*((distance_matrix)./theta)).*exp(-sqrt(3).*((distance_matrix)./theta));
    case 'Matern_5_2'
        K = @(theta, distance_matrix) (1 + sqrt(5).*(distance_matrix./theta) ...
            + (5/3).*((distance_matrix./theta).^2)).*exp(-sqrt(5)*(distance_matrix./theta));
    otherwise
        error('I cannot find this correlation family.')
end

dim = length(theta0);
R_aux = ones(size(scaled_X, 1), size(scaled_X, 1));

parfor i = 1:dim
    distance_matrix = squareform(pdist(scaled_X(:, i)));
    R = K(theta(i), distance_matrix);
    R_aux = R.*R_aux;
end

R = R_aux;
R = R + (10+length(R))*1e-10*eye(length(R));
[C, flag] = chol(R, 'lower');
F = ones(length(R),1);
if flag ~= 0, error('Correlation matrix is NOT symmetric positive definite'); end
metamodel_pars.fit.F_tilde = C\F; %inv(C)*F
metamodel_pars.fit.Y_tilde = C\(scaled_Y);
[Q, G ] = qr(metamodel_pars.fit.F_tilde,0);

if rcond(G) <1e-10
    if cond(F) > 1e-15
        error('F is too ill conditioned. Poor combination of regression model and design sites');
    end
end

% Return fitting parameters
metamodel_pars.fit.beta = G \ (Q'*metamodel_pars.fit.Y_tilde);
metamodel_pars.fit.sigmaSQ = (1/length(R))*(norm(metamodel_pars.fit.Y_tilde - metamodel_pars.fit.F_tilde*metamodel_pars.fit.beta))^2;
metamodel_pars.fit.gama = transpose(inv(C))*(metamodel_pars.fit.Y_tilde - metamodel_pars.fit.F_tilde*metamodel_pars.fit.beta);
metamodel_pars.fit.invC = inv(C);
metamodel_pars.fit.invG = inv(G);
metamodel_pars.fit.phi = prod(diag(C)).^(2/length(R))*metamodel_pars.fit.sigmaSQ;
metamodel_pars.fit.family = metamodel_pars.correlation_families{idx_family};
end

