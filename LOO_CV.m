function [ all_families ] = LOO_CV( all_families )

% [ shuffe_position, all_families] = shuffle_data( all_families );

for i = 1:length(all_families)
    
    aux_X = all_families(i).scaled_X;
    aux_Y = all_families(i).scaled_Y;
    
    for j = 1:length(aux_Y)
        metamodel_pars = all_families(i);
%         metamodel_pars.Input.X_coordinates(j,:) = [];
%         metamodel_pars.Input.Y_coordinates(j,:) = [];
        metamodel_pars.scaled_X(j,:) = [];
        metamodel_pars.scaled_Y(j,:) = [];
        metamodel_pars = objective_function( metamodel_pars, i);
        metamodel_pars = asses_r0(metamodel_pars, num2cell(aux_X(j,:)));
        aux_Y(j,2) = metamodel_pars.fit.beta + metamodel_pars.r0(1,:)*metamodel_pars.fit.gama;
    end
    all_families(i).MSE = (1/length(aux_X)).*sum((aux_Y(:,1) - aux_Y(:,2)).^2);
    all_families(i).RMSE = sqrt(all_families(i).MSE);
    all_families(i).NRMSE = all_families(i).RMSE/(max(aux_Y(:,1))-min(aux_Y(:,2)));
end

% [ all_families] = organize_data( all_families, shuffe_position);

end

