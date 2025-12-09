function FIGURE7_localSensitivityAnalysis(regenerate)
% This function performs local sensitivity analysis around the parameter
% values obtained via running the function RUNNER_DunABC. That is, it takes
% the particle set of posterior samples and evaluates the parameter
% sensitivity (of the steady state) at each
%
% Steady state sensitivity is defined by considering the following:
%
% For an ODE model x'(t) = F(x; theta), steady states must satisfy
%
%                   F(x; theta)   = 0
%         d/dtheta( F(x; theta) ) = 0
%     dF/dx dx/dtheta + dF/dtheta = 0
%                       dx/dtheta = - (dF/dx)^-1 dF/dtheta
%                                 = - (Jx)^-1 Jtheta           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify whether to calculate relative (%) sensitivities using logs, or
% raw sensitivities that are affected by scale
use_log_sensitivity = true;

% Define the names of the variables and the parameters
variable_names = {'FS', 'fs', 'SL', 'sl', 'I'};
parameter_names = {'$\lambda_{F\!S}$', '$K_{F\!S}$', '$\lambda_{f\!s}$', '$\lambda_{S\!L}$', '$\lambda_{sl}$', '$\lambda_{I}$', '$\alpha_{d}$', '$\alpha_{u}$'};

% For now, sensitivity is considered for the wild-type parameterisation
phi_s = [1;1;1];
phi_r = [1;1];

% Filename specification
filename = 'sensitivity_data.mat';

% Specifies what percentage of results to remove in plotting credible
% intervals (e.g. 0.05 gives the 2.5% - 97.5% interval)
plot_quantile = 0.05;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the particles
load('final_particles.mat','particles','prior');
% Read out their discrepancies
part_Ds = getProperty(particles,'D');


% Check if results are to be regenerated
if nargin < 1
    regenerate = false;
end
if regenerate || ~exist(filename, 'file')
    
    %%% SENSITIVITY ACROSS DIFFERENT PARTICLES
    
    % Count how many particles were provided
    Nparts = length(particles);
    
    % Just use a hard-coded specfication of problem dimensions for now
    Nx = 5;
    Ntheta = 8;
    
    % Initialise sensitivity matrix array (one per particle)
    S = nan(Nx, Ntheta, Nparts);
    % Ditto for prior samples
    Sprior = nan(Nx, Ntheta, Nparts);
    
    % Iterate over particles, processing the local sensitivity analysis for
    % each in turn
    for p = 1:Nparts
        
        % Read out the current particle's parameters, including (hard-coded)
        % conversion of log-rates to rates
        params = particles{p}.theta;
        params(1:6) = 10.^(params(1:6));
        
        % Run the model to find network's steady state for these parameters
        [Xstar, converged] = DunSS_paramODE(params, phi_s, phi_r);
        
        % The calibrated particles should converge, but check for safety
        if converged
            
            % Calculate the Jacobians, dF/dx and dF/dtheta
            [Jvar, Jparam] = calcJacobians(Xstar, params, phi_s, phi_r);
            
            % Determine sensitivity of the steady state to the parameters
            % for this particle and store in the sensitivity array
            if use_log_sensitivity
                % NOTE: Here pre-multiplication by diag(1/Xstar) and
                % post-multiplication by diag(params) is used to make the
                % sensitviity values calculated percentage change in steady
                % state value for a percentage change in parameter value
                S(:,:,p) = -diag(1./Xstar(1:Nx)) * (Jvar \ Jparam) * diag(params);
            else
                S(:,:,p) = -Jvar \ Jparam;
            end
            
        % If we do hit the failure case, output a warning to the user and 
        % leave the sensitivity matrix as NaN's for this particle
        else
            warning('A particle failed to converge to steady state when simulated!');
        end
        
    end
    
    % Iterate over prior samples as well, processing the local sensitivity
    % analysis for each in turn
    for p = 1:Nparts
        
        % Generate a sample from the prior (uniform)
        prior_sample = prior.sample();
        % Exponentiate the rate parameters and half-saturation constant (log-uniform prior on these)
        prior_sample(1:6) = 10.^(prior_sample(1:6));
        
        % Run the model to find network's steady state for these parameters
        [Xstar, converged] = DunSS_paramODE(prior_sample, phi_s, phi_r);
        
        % The calibrated particles should converge, but check for safety
        if converged
            
            % Calculate the Jacobians, dF/dx and dF/dtheta
            [Jvar, Jparam] = calcJacobians(Xstar, prior_sample, phi_s, phi_r);
            
            % Determine sensitivity of the steady state to the parameters
            % for this particle and store in the sensitivity array
            if use_log_sensitivity
                % NOTE: Here pre-multiplication by diag(1/Xstar) and
                % post-multiplication by diag(params) is used to make the
                % sensitviity values calculated percentage change in steady
                % state value for a percentage change in parameter value
                Sprior(:,:,p) = -diag(1./Xstar(1:Nx)) * (Jvar \ Jparam) * diag(prior_sample);
            else
                Sprior(:,:,p) = -Jvar \ Jparam;
            end
            
        % If we do hit the failure case, output a warning to the user and 
        % leave the sensitivity matrix as NaN's for this particle
        else
            warning('A particle failed to converge to steady state when simulated!');
        end
        
    end
    
    % Save the sensitivity data
    save(filename, 'S', 'Sprior');
    
else
    load(filename, 'S', 'Sprior');
end

% Use the sensitivity matrix to establish dimensions
[Nx, Ntheta, Nparts] = size(S);
   
% Calculate statistical properties across set of sensitivity matrices found
Smean = mean(S,3,'omitnan');
Sstd = std(S,[],3,'omitnan');



%%% SENSITIVITY OF THE BASE MODEL

% Base parameters - unity for everything, except 0.5 for transport terms
pbase = [1,1,1,1,1,1,0.5,0.5];

% Run the base model and warn if it didn't converge
[Xbase, converged] = DunSS_paramODE(pbase, phi_s, phi_r);
if ~converged
    warning('The base model did not converge! Comparisons to base model sensitivity are likely incorrect.');
end

% Calculate the sensivitiy using the Jacobians, as above
[Jvar, Jparam] = calcJacobians(Xbase, pbase, phi_s, phi_r);
if use_log_sensitivity
    Sbase = -diag(1./Xbase(1:Nx)) * (Jvar \ Jparam) * diag(pbase);
else
    Sbase = -Jvar \ Jparam; 
end

% Compare the number of sign mismatches between the mean and median
% sensitivity values, and the base parameter values
fprintf('Percentage of sensitivity directions in agreement between mean sensitivity and base sensitivity is %g%%\n', 100 - mean(abs(sign(Sbase) - sign(Smean)),'all')*100);

% Also count out what proportion of the total particle set show all signs
% in agreement
agree_count = 0;
for k = 1:Nparts
    
    % Grab out sensitivity for this particle
    Spart = squeeze(S(:,:,k));
    % Increment counter if all signs in agreement
    if sum(abs(sign(Sbase) - sign(Spart)),'all') == 0
        agree_count = agree_count + 1;
    end
    
end
fprintf('The Percentage of particles with all signs in agreement is: %g%%\n', agree_count / Nparts * 100);
    
% Calculate and output the "Z-score" matrix
Z = (Smean - Sbase) ./ Sstd



%%% SENSITIVITY OF THE BEST MODEL

% Best parameters - grab out minimum discrepancy particle
[~,minloc] = min(part_Ds);
pbest = particles{minloc}.theta;
pbest = [exp( pbest(1:6) ), pbest(7:8)]

% Run the base model and warn if it didn't converge
[Xbest, converged] = DunSS_paramODE(pbest, phi_s, phi_r);
if ~converged
    warning('The best model did not converge!');
end

% Calculate the sensivitiy using the Jacobians, as above
[Jvar, Jparam] = calcJacobians(Xbest, pbest, phi_s, phi_r);
if use_log_sensitivity
    Sbest = -diag(1./Xbest(1:Nx)) * (Jvar \ Jparam) * diag(pbest);
else
    Sbest = -Jvar \ Jparam; 
end

% Compare the number of sign mismatches between the mean and median
% sensitivity values, and the base parameter values
fprintf('Percentage of sensitivity directions in agreement between best sensitivity and base sensitivity is %g%%\n', 100 - mean(abs(sign(Sbest) - sign(Sbase)),'all')*100);
 


%%% SENSITIVITY OF ONE OF THE "EXTREME SENSITIVITY" MODELS

% Pick out a random example (first one satisfying a simple condition that
% identifies these models) of a calibrated model with extreme sensitivity
plus_loc = find( squeeze(S(1,3,:)) < -1.9, 1);
Splus = squeeze(S(:,:,plus_loc));

% Also output the parameters to the user (for writing numerical parameter
% values in the document)
pplus = particles{plus_loc}.theta;
pplus = [ 10.^pplus(1:6), pplus(7:8)]




%%% GRAPHICAL DISPLAYS OF SENSITIVITY

% Prepare a red to blue colormap for negative/positive sensitivity
blue_red = [ linspace(0.6,1.0,50)', linspace(0.0,1.0,50)', linspace(0.0,1.0,50)';
             linspace(1.0,0.0,50)', linspace(1.0,0.0,50)', linspace(1.0,0.6,50)'  ];

% Determine the maximum sensitivity value across both matrices
max_sens = max([Smean(:); Sbase(:); Sbest(:)]);


%%% Mean sensitivity

% Prepare figure
figure('units','Normalized','OuterPosition',[0.05 0.2 0.5 0.5]);
hold on;

% Calculate the relative sensitivities of each parameter (columns)
param_sensitivities = sum( abs(Smean), 1 );
param_sensitivities = param_sensitivities / max(param_sensitivities);
% Calculate the relative sensitivities of each parameter (columns)
X_sensitivities = sum( abs(Smean), 2 );
X_sensitivities = X_sensitivities / max(X_sensitivities);

% Convert these into cell arrays of strings
X_sens_txts = cell(1,Nx);
for k = 1:Nx
    X_sens_txts{k} = num2str(X_sensitivities(k),'%.2f');
end
param_sens_txts = cell(1,Ntheta);
for k = 1:Ntheta
    param_sens_txts{k} = num2str(param_sensitivities(k),'%.2f');
end

% Plot the matrix
% (Add a column and row to the matrix to deal with pcolor shenanigans)
plot_mat = [Smean, zeros(Nx,1); zeros(1,Ntheta), 0];
pcolor(plot_mat);
caxis(max_sens*[-1 1]);
colormap(blue_red);

% Loop over each element, adding labels with numerical value
for i = 1:size(Smean,1)
    for j = 1:size(Smean,2)
        % Use black text for weak sensitivity (light background)
        if abs(Smean(i,j)) < max_sens/2
            text(j+0.5,i+0.5, num2str(Smean(i,j),'%.2f'), 'Color', [0 0 0], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        % Use white text for strong sensitivity (dark background)
        else
            text(j+0.5,i+0.5, num2str(Smean(i,j),'%.2f'), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        end
    end
end

% Fix axis limits
axis equal;
xlim([1 Ntheta+1]);
ylim([1 Nx+1]);

% Adjust axis properties and add labels
set(gca,'FontSize',20,'XAxisLocation','top','YDir','reverse','TickLabelInterpreter','LaTeX');
ax_obj = gca;
ax_obj.XRuler.TickLabelGapOffset = -5;
ax_obj.YRuler.TickLabelGapOffset = -1.5;
xticks(0.5+(1:Ntheta));
yticks(0.5+(1:Nx));
xticklabels(parameter_names);
yticklabels(variable_names);
% Add title
title('Parameterised (Mean)','FontSize',22);


%%% Base sensitivity

% Prepare figure
figure('units','Normalized','OuterPosition',[0.05 0.2 0.5 0.5]);
hold on;

% Calculate the relative sensitivities of each parameter (columns)
param_sensitivities = sum( abs(Sbase), 1 );
param_sensitivities = param_sensitivities / max(param_sensitivities);
% Calculate the relative sensitivities of each parameter (columns)
X_sensitivities = sum( abs(Sbase), 2 );
X_sensitivities = X_sensitivities / max(X_sensitivities);

% Convert these into cell arrays of strings
X_sens_txts = cell(1,Nx);
for k = 1:Nx
    X_sens_txts{k} = num2str(X_sensitivities(k),'%.2f');
end
param_sens_txts = cell(1,Ntheta);
for k = 1:Ntheta
    param_sens_txts{k} = num2str(param_sensitivities(k),'%.2f');
end

% Plot the matrix
% (Add a column and row to the matrix to deal with pcolor shenanigans)
plot_mat = [Sbase, zeros(Nx,1); zeros(1,Ntheta), 0];
pcolor(plot_mat);
caxis(max_sens*[-1 1]);
colormap(blue_red);

% Loop over each element, adding labels with numerical value
for i = 1:size(Sbase,1)
    for j = 1:size(Sbase,2)
        % Use black text for weak sensitivity (light background)
        if abs(Sbase(i,j)) < max_sens/2
            text(j+0.5,i+0.5, num2str(Sbase(i,j),'%.2f'), 'Color', [0 0 0], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        % Use white text for strong sensitivity (dark background)
        else
            text(j+0.5,i+0.5, num2str(Sbase(i,j),'%.2f'), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        end
    end
end

% Fix axis limits
axis equal;
xlim([1 Ntheta+1]);
ylim([1 Nx+1]);

% Adjust axis properties and add labels
set(gca,'FontSize',20,'XAxisLocation','top','YDir','reverse','TickLabelInterpreter','LaTeX');
ax_obj = gca;
ax_obj.XRuler.TickLabelGapOffset = -5;
ax_obj.YRuler.TickLabelGapOffset = -1.5;
xticks(0.5+(1:Ntheta));
yticks(0.5+(1:Nx));
xticklabels(parameter_names);
yticklabels(variable_names);
% Add title
title('Parameter-Free','FontSize',22);


%%% Best sensitivity

% Prepare figure
figure('units','Normalized','OuterPosition',[0.05 0.2 0.5 0.5]);
hold on;

% Calculate the relative sensitivities of each parameter (columns)
param_sensitivities = sum( abs(Sbest), 1 );
param_sensitivities = param_sensitivities / max(param_sensitivities);
% Calculate the relative sensitivities of each parameter (columns)
X_sensitivities = sum( abs(Sbest), 2 );
X_sensitivities = X_sensitivities / max(X_sensitivities);

% Convert these into cell arrays of strings
X_sens_txts = cell(1,Nx);
for k = 1:Nx
    X_sens_txts{k} = num2str(X_sensitivities(k),'%.2f');
end
param_sens_txts = cell(1,Ntheta);
for k = 1:Ntheta
    param_sens_txts{k} = num2str(param_sensitivities(k),'%.2f');
end

% Plot the matrix
% (Add a column and row to the matrix to deal with pcolor shenanigans)
plot_mat = [Sbest, zeros(Nx,1); zeros(1,Ntheta), 0];
pcolor(plot_mat);
caxis(max_sens*[-1 1]);
colormap(blue_red);

% Loop over each element, adding labels with numerical value
for i = 1:size(Sbest,1)
    for j = 1:size(Sbest,2)
        % Use black text for weak sensitivity (light background)
        if abs(Sbest(i,j)) < max_sens/2
            text(j+0.5,i+0.5, num2str(Sbest(i,j),'%.2f'), 'Color', [0 0 0], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        % Use white text for strong sensitivity (dark background)
        else
            text(j+0.5,i+0.5, num2str(Sbest(i,j),'%.2f'), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        end
    end
end

% Fix axis limits
axis equal;
xlim([1 Ntheta+1]);
ylim([1 Nx+1]);

% Adjust axis properties and add labels
set(gca,'FontSize',20,'XAxisLocation','top','YDir','reverse','TickLabelInterpreter','LaTeX');
ax_obj = gca;
ax_obj.XRuler.TickLabelGapOffset = -5;
ax_obj.YRuler.TickLabelGapOffset = -1.5;
xticks(0.5+(1:Ntheta));
yticks(0.5+(1:Nx));
xticklabels(parameter_names);
yticklabels(variable_names);
% Add title
title('Parameterised (Best)','FontSize',22);




%%% Supplement - extreme sensitivity of many individual calibrations

% Prepare figure
figure('units','Normalized','OuterPosition',[0.05 0.2 0.5 0.5]);
hold on;

% Plot the matrix
% (Add a column and row to the matrix to deal with pcolor shenanigans)
plot_mat = [Splus, zeros(Nx,1); zeros(1,Ntheta), 0];
pcolor(plot_mat);
% Rescale maximum sensitivity for this plot
max_sens = max(abs(Splus(:)));
caxis(max_sens*[-1 1]);
colormap(blue_red);

% Loop over each element, adding labels with numerical value
for i = 1:size(Splus,1)
    for j = 1:size(Splus,2)
        % Use black text for weak sensitivity (light background)
        if abs(Splus(i,j)) < max_sens/2.1
            text(j+0.5,i+0.5, num2str(Splus(i,j),'%.2f'), 'Color', [0 0 0], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        % Use white text for strong sensitivity (dark background)
        else
            text(j+0.5,i+0.5, num2str(Splus(i,j),'%.2f'), 'Color', [1 1 1], 'HorizontalAlignment', 'Center', 'VerticalAlignment','Middle', 'FontSize', 16)
        end
    end
end

% Fix axis limits
axis equal;
xlim([1 Ntheta+1]);
ylim([1 Nx+1]);

% Adjust axis properties and add labels
set(gca,'FontSize',20,'XAxisLocation','top','YDir','reverse','TickLabelInterpreter','LaTeX');
ax_obj = gca;
ax_obj.XRuler.TickLabelGapOffset = -5;
ax_obj.YRuler.TickLabelGapOffset = -1.5;
xticks(0.5+(1:Ntheta));
yticks(0.5+(1:Nx));
xticklabels(parameter_names);
yticklabels(variable_names);
% Add title
title('Extreme Sensitivity Model','FontSize',22);




%%% Mean Parameter Sensitivity Plot

% Calculate sensitivity (average absolute value) for each parameter across
% all calibrated models
all_param_sens = squeeze(mean(abs(S), 1))';
prior_sens = squeeze(mean(abs(Sprior), 1))';

% Extract 95% credible intervals for these
all_param_min_sens = quantile(all_param_sens,plot_quantile/2);
all_param_max_sens = quantile(all_param_sens,1 - plot_quantile/2);
prior_min_sens = quantile(prior_sens,plot_quantile/2);
prior_max_sens = quantile(prior_sens,1 - plot_quantile/2);

% Calculate sensitivity for param-free and best models, too
base_param_sens = mean(abs(Sbase));
best_param_sens = mean(abs(Sbest));

% Create also a plot of overall sensitivity to the parameters
figure('units','Normalized','OuterPosition',[0.05 0.2 0.4 0.5]);
hold on;
for p = 1:Ntheta
    
    plot([p p], [prior_min_sens(p), prior_max_sens(p)], 'LineWidth',5, 'Color', [0.65 0.65 0.65]);
    plot([p-0.2 p+0.2], [prior_min_sens(p), prior_min_sens(p)], 'LineWidth',5, 'Color', [0.65 0.65 0.65]);
    plot([p-0.2 p+0.2], [prior_max_sens(p), prior_max_sens(p)], 'LineWidth',5, 'Color', [0.65 0.65 0.65]);
    plot([p p], [all_param_min_sens(p), all_param_max_sens(p)], 'k', 'LineWidth',5);
    plot([p-0.2 p+0.2], [all_param_min_sens(p), all_param_min_sens(p)], 'k', 'LineWidth',5);
    plot([p-0.2 p+0.2], [all_param_max_sens(p), all_param_max_sens(p)], 'k', 'LineWidth',5);
    plot([p-0.2 p+0.2], [base_param_sens(p), base_param_sens(p)], 'r', 'LineWidth',5);
    plot([p-0.2 p+0.2], [best_param_sens(p), best_param_sens(p)], 'b', 'LineWidth',5);
    
    % Add legend using second plotted data (simply to put credible interval
    % last in legend, but first in plot)
    if p == 2
        leg_obj = legend({'','','','','','','Parameter-free', 'Best Parameterised', '95% CI (Prior)', '', '', '95% CI (Calibrated)'},'FontSize',20);
        leg_obj.AutoUpdate = false;
    end
        
end

xlim([0.5, Ntheta+0.5]);
xticklabels(parameter_names);
ax_obj = gca;
ax_obj.FontSize = 22;
ax_obj.XLabel.FontSize = 24;
ax_obj.TickLabelInterpreter = 'LaTeX';

end