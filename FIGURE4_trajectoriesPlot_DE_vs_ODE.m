function FIGURE4_trajectoriesPlot_DE_vs_ODE
% This function generates Figure 4 in the paper, which compares the
% trajectories taken by the difference equation formulation of the model to
% those taken by the ODE version of the model

%%% UNIVERSAL SPECIFICATIONS

% Define the initial condition
X0 = 1.15*ones(6,1);

% Specify whether to use the delayed model for the difference equation
use_delayed = false;

% Values to use for alpha_down and alpha_up
alphas = [0.7, 0.75];

% Figure outer position (window size)
figOuterPosition = [0.09, 0.12, 0.72, 0.76];


%%% FIRST PANEL: A RUN OF THE WILD TYPE PLANT

% Define the plotting colours
plot_colors = [ [0.4, 0.4, 1.0];
                [0.0, 0.0, 0.6];
                [1.0, 0.4, 0.4];
                [0.6, 0.0, 0.0];
                [0.8, 0.0, 0.8];
                [0.5, 0.5, 0.5] ];

% Number of iterations to run for the first scenario
N_iters = 20;

% Set genotypes to the wild type
phi_s = [1;1;1];
phi_r = [1;1];

% Number of legend columns
num_cols = 3;



% Store the initial condition in the state variable vector (and blank
% values to fill in after). Add values for FS_old and sl_old if delayed
if use_delayed
    X = [X0;X0(1);X0(4)];
else
    X = X0;
end
X_DE = [X, nan(size(X,1),N_iters)];

% Run the model in its difference equation form
for i = 1:N_iters
    
    % Update the value of the current state
    X = DunDifference(X, alphas, phi_s, phi_r, use_delayed);
    
    % Store the updated state
    X_DE(:,i+1) = X;
    
end

% Run the model in its ODE form
[T, X_ODE] = ode15s(@(t,X) DunODE(t, X, alphas, phi_s, phi_r), [0 N_iters], X0);
X_ODE = X_ODE';

% Prepare figure
figure('units','normalized','OuterPosition',figOuterPosition);
hold on;

% Plot the two against one another
for k = 1:size(X_ODE,1)
    
    % Plot the ODE as a solid line
    plot(T, X_ODE(k,:), 'LineWidth',3, 'Color', plot_colors(k,:) );
    
    % Plot the difference equation as a dotted line between marked points
    % at each iteration's value
    plot(0:N_iters, X_DE(k,:), '.:', 'MarkerSize',40,'LineWidth',1.5,'Color',plot_colors(k,:));
    
end

% Axis clean-up
axis square;
xlabel('Time/Number of Iterations');
ylabel('Dimensionless Concentration');
set(gca, 'FontSize',20, 'LineWidth',2.5);
    
% Add legend
species_names = {'[FS]', '[fs]', '[SL]', '[sl]', '[I]', '[ck]'};
leg_txt = {};
for k = 1:length(species_names)
    % Add entry (with spaces if not final set)
    if k < length(species_names) * (num_cols / (num_cols + 1)) + 1e-8
        leg_txt{end+1} = [' ',species_names{k},'   {\color{white}x}'];
    else
        leg_txt{end+1} = [' ',species_names{k}];
    end
    % Add blank entry because we have two plots per
    leg_txt{end+1} = '';
end
leg_obj = legend(leg_txt);
leg_obj.FontSize = 20;
leg_obj.NumColumns = num_cols;
leg_obj.Box = 'off';


%%% SECOND PANEL: DIFFERENT GENOTYPES

% Number of iterations to use for this panel
N_iters = 60;

% Define the plotting colours
plot_colors = [ [0.0, 0.0, 0.0];
                [0.0, 0.6, 0.0];
                [0.3, 0.8, 0.3];
                [0.0, 0.9, 0.7] ];

% Use the special value for alphas that causes instability in one case
alphas = [1,1];
            
% Listing genotypes to plot - wild type
genotype_list{1}.phi_s = [1;1;1];
genotype_list{1}.phi_r = [1;1];
% Listing genotypes to plot - no root-derived feedback signal
genotype_list{2}.phi_s = [1;1;1];
genotype_list{2}.phi_r = [0;1];
% Listing genotypes to plot - no root-derived strigolactone
genotype_list{3}.phi_s = [1;1;1];
genotype_list{3}.phi_r = [1;0];
% Listing genotypes to plot - no shoot-derived strigolactone
%      (special case for which this is unstable)
genotype_list{4}.phi_s = [1;0;1];
genotype_list{4}.phi_r = [0;1];

% Number of legend columns
num_cols = 2;
    
% Prepare figure
figure('units','Normalized','OuterPosition',figOuterPosition);
hold on;

% Loop over each genotype in the list
for g = 1:length(genotype_list)
    
    % Read out genotype 
    phi_s = genotype_list{g}.phi_s;
    phi_r = genotype_list{g}.phi_r;
    
    
    %%% Simulate the difference equation model
    
    % Prepare initial condition
    if use_delayed
        X = [X0;X0(1);X0(4)];
    else
        X = X0;
    end
    X_DE = [X, nan(size(X,1),N_iters)];

    % Iterate model right hand side and store
    for i = 1:N_iters
        X = DunDifference(X, alphas, phi_s, phi_r, use_delayed);
        X_DE(:,i+1) = X;
    end
    
    
    %%% Simulate the ODE model
    
    [T, X_ODE] = ode15s(@(t,X) DunODE(t, X, alphas, phi_s, phi_r), [0 N_iters], X0);
    X_ODE = X_ODE';
    
    
    %%% Add the results to the plot
    
    % Plot the shoot feedback signal for both model frameworks
    plot(T, X_ODE(1,:), 'LineWidth', 3, 'Color', plot_colors(g,:) );
    plot(0:N_iters, X_DE(1,:), '.:', 'MarkerSize', 40, 'LineWidth', 1.5, 'Color', plot_colors(g,:) );
    
    
end
    
% Axis clean-up
axis square;
xlabel('Time/Number of Iterations');
ylabel('Concentration, [FS] (Dimensionless)');
set(gca, 'FontSize',20, 'LineWidth',2.5);

% Add legend
genotype_names = {'Wild Type', 'No Root FS', 'No Root SL', 'No Shoot SL'};
leg_txt = {};
for k = 1:length(genotype_names)
    % Add entry (with spaces if not final set)
    if k < length(genotype_names) * (num_cols / (num_cols+1)) + 1e-6
        leg_txt{end+1} = [' ',genotype_names{k},'   {\color{white}x}'];
    else
        leg_txt{end+1} = [' ',genotype_names{k}];
    end
    % Add blank entry because we have two plots per
    leg_txt{end+1} = '';
end
leg_obj = legend(leg_txt);
leg_obj.FontSize = 20;
leg_obj.NumColumns = 2;
leg_obj.Box = 'off';