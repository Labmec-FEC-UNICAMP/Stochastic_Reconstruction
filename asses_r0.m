function [ metamodel_pars ] = asses_r0( metamodel_pars, scaled_post_processing_grid )


switch metamodel_pars.fit.family
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

dim = length(metamodel_pars.theta0);

if iscell(scaled_post_processing_grid)
    idx = cell(1, dim);
    for k = 1:dim
        idx{k} = reshape(scaled_post_processing_grid{k}, [numel(scaled_post_processing_grid{k}) 1]);
    end
    plan_grid = cell2mat(idx);
    
else
    plan_grid = scaled_post_processing_grid;
end

r_aux = zeros(size(plan_grid, 1), size(metamodel_pars.scaled_X,1));
r = zeros(size(metamodel_pars.scaled_X, 1), dim);


for j = 1:size(plan_grid, 1)
    for i = 1:dim
        distance_matrix = abs(metamodel_pars.scaled_X(:,i) - plan_grid(j,i));
        r(:,i) = K(metamodel_pars.theta(i), distance_matrix);
        %         r_aux(j,:) = transpose(r.*r_aux(j,:)');
    end
    r_aux(j, :) = transpose(prod(r,2));
end

metamodel_pars.r0 = r_aux;
metamodel_pars.plan_grid = plan_grid;
end

