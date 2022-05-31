fclose('all');
close all;
clear;
clc;


%% Input Data

%get data
% columns are organized as follows:
% coord_x | coord_y | strain_xx | strain_yy | strain_xy

data = 'data_Sample21.csv';
data_read = csvread(data, 0, 0);

%%

A = [data_read(:,1) data_read(:,2)]; %coordinates
AUX = [data_read(:,3) data_read(:,4) data_read(:,5)]; %strains

% Test each strain
for i = 1:3
    
    % Organize and normalize data *check with Anderson!
    A = [A, AUX(:,i)];
    X = A(:, 1:(end-1));
    Y = A(:, end);
    Y = (Y - mean(Y))./std(Y); %Transforma??o iso-probabilistica
    
    % X = metamodel_pars.Input.X_coordinates;
    % Y = metamodel_pars.Input.Y_coordinates;
    metamodel_pars.scaled_X = X;
    metamodel_pars.scaled_Y = Y;
    
    
    %%
    % Firt guess
    metamodel_pars.theta0 = 0.05*ones(1, size(metamodel_pars.scaled_X,2));
    metamodel_pars.theta_lower = 0.01*ones(1, size(metamodel_pars.scaled_X,2));
    metamodel_pars.theta_upper = 1*ones(1, size(metamodel_pars.scaled_X,2));
    
    % Optimization algorithm by Lophaven 2002
    [~, all_families] = determine_correlation_family(metamodel_pars);
    %%
    
    [ all_families ] = LOO_CV( all_families );
    [~, aux ] = min([all_families(:).MSE]);
    metamodel_pars = all_families(aux);
    %%
    
    if i == 1
        fprintf('Strain_xx kernel is computed')
        resultados_xx = all_families;
    elseif i == 2
        fprintf('Strain_yy kernel is computed')
        resultados_yy= all_families;
    else
        fprintf('Strain_xy kernel is computed')
        resultados_xy = all_families;
    end
    A(:,end) = [];
end

