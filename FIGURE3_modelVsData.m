function FIGURE3_modelVsData
% This function produces the figure presenting comparisons between the best
% particle's model, the baseline (parameter-free) model, and the data from
% experiments. Given the variability in scale of model predictions, these
% predictions are mapped onto values similar to the discrete values
% recorded for the experiments using a linear transformation that preserves
% rank and thus does not affect the discrepancy measure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify whether to label discrepancies
label_discrepancies = true;

% Specify maximum number of experiments to put on a single plot
plot_max = 26;

% Specify bar colours for the plots
bar_colors = [0.0, 0.0, 0.0;        % Experiments - black
              1.0, 0.2, 0.2;        % Baseline - red
              0.2, 0.2, 1.0  ];     % Calibrated - blue

% List the species contained in the experiments
species = {'Cytokinin','RMS1-shoot','RMS1-root','Branching Inhibition'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% INITIAL PREPARATION

% Load in the particle data
load('final_particles.mat','particles');
% Load in the experimental data
experiments = DunExperiments;

% Pull out the lowest-discrepancy particle's results
part_Ds = getProperty(particles,'D');
[~,loc] = min(part_Ds);
model_y = particles{loc}.y;

% Also evaluate the parameter-free model
default_params = [1, 1, 1, 1, 1, 1, 0.5, 0.5];

% Generate predictions for the default model
default_y = runAllExperiments(default_params, experiments);

% Store the two models' predictions together
predictions = {default_y, model_y};



%%% PLOTTING FUNCTION

% Loop over the different experiments, processing each separately
for k = 1:length(experiments)
    
    % Read out the experimental results for this experiment
    exp_y = experiments{k}.results;
    
    % Sort results so that all two-component experiments are at the end
    [~,I] = sort(~isnan(exp_y(:,2)));
    exp_y = exp_y(I,:);
    % Set value for start of two-component experiments (-1 if not present)
    twocomp_start = -1 + ~isempty(find(~isnan(exp_y(:,2)),1)) * (size(exp_y,1)+1);
    
    % Now grab only the actual observations and compile into column vector
    template = ~isnan(exp_y);
    exp_y = exp_y(template);
    
    % Apply re-ordering and template to the descriptions and genotype info
    exp_mutants = experiments{k}.mutants(I,:);
    exp_mutants = exp_mutants(template)';
    
    % Also grab out genotype information in re-arranged order
    exp_phi_s = experiments{k}.phi_s(I);
    exp_phi_r = experiments{k}.phi_r(I);
    % Ensure these are column cell arrays, and duplicate them
    exp_phi_s = [exp_phi_s(:),exp_phi_s(:)];
    exp_phi_r = [exp_phi_r(:),exp_phi_r(:)];
    % Apply the template
    exp_phi_s = exp_phi_s(template);
    exp_phi_r = exp_phi_r(template);
        
    % Loop over each set of model predictions, and scale each in turn
    model_fits = zeros(length(exp_y), length(predictions));
    for p = 1:length(predictions)
        
        % Read out the model's predictions for this experiment
        model_y = predictions{p}{k};
        
        % Apply same re-ordering applied to data to the model output
        model_y = model_y(I,:);
        
        % Grab out only the valid ones using the template
        model_y = model_y(template);
        
        % Select the optimal scaling for these values to match the data.
        % Linear models with a variety of link functions are trialled
        model_y_trans = robustLinearTransform(model_y, exp_y);
        
        % Store this in the bar graph structure
        model_fits(:,p) = model_y_trans;
        
    end
    
    % Find discrepancies between experimental results and model predictions
    for p = 1:length(predictions)
        discrep_flags{p} = findMismatches(model_fits(:,p),exp_y);
    end
    
    % Space-saving measures for the case with many measurements
    if strcmp(species{k},'Branching Inhibition')
        
        % Cut down labels so they occupy less space
        exp_mutants = replace(exp_mutants,'rms','r');
        exp_mutants = replace(exp_mutants,'WT','W');
        
        % Identify 'null' experiments guaranteed to not produce any
        % branching inhibition (knockout of RMS3/RMS4 in shoot)
        c2 = 0;
        null = [];
        for e = 1:length(exp_phi_s)
            
            % Check if this plant produces no branch inhibition (shoot)
            % or no strigolactones at all (shoot and root)
            if all(exp_phi_s{e}(3,:) == 0) || ( all(exp_phi_s{e}(2,:) == 0) && all(exp_phi_r{e}(2,:) == 0) )
                % Add this value to the deletion list
                c2 = c2 + 1;
                null(c2) = e;
            end
            
        end
        
        % Remove these null experiments
        exp_y(null) = [];
        exp_mutants(null) = [];
        model_fits(null,:) = [];
        for p = 1:length(predictions)
            discrep_flags{p}(null) = [];
        end
        
        % Shift the start of two-component results based on deletions
        twocomp_start = twocomp_start - sum(null < twocomp_start);
        
    end
    
    % Plot the bar graph(s) for this experiment
    plotting_done = false;
    c = 0;
    while ~plotting_done
        
        % Increment counter
        c = c + 1;
        % Increment index
        ind = (1:plot_max) + plot_max*(c-1);
        % When hitting end of data, trim to length of data and stop looping
        if max(ind) >= length(exp_y)
            ind(ind > length(exp_y)) = [];
            plotting_done = true;
        end
        
        % Italicise all mutant names
        for r = 1:5
            exp_mutants = replace(exp_mutants,['rms',num2str(r)],['{\itrms',num2str(r),'}']);
            exp_mutants = replace(exp_mutants,['r',num2str(r)],['{\itr',num2str(r),'}']);
        end
        
        % Create a new figure
        figure('units','Normalized','OuterPosition',[0 0 1 1]);
        
        % Plot the bar graph for this portion of the data
        bar_obj = bar([exp_y(ind),model_fits(ind,:)]);
        % Add a tick for each plotted bar
        xticks(1:length(ind));
        % Label each bar with the mutant
        xticklabels(exp_mutants(ind));
        % Set the y-axis limit
        ylim([0 max(exp_y)*1.1]);
        % Add ticks if this is the first plot for this set, otherwise blank
        if c == 1
            yticks(1:max(exp_y));
            % Also label the measured quantity on first plot of set
            ylabel(species{k},'FontSize',20);
        else
            set(gca,'ytick',[]);
            set(gca,'yticklabels',[]);
        end
        
        % If this is not the first plot, also delete the right axis so
        % plots can be joined together later
        if c>1
            ax_obj = gca;
            ax_obj.YAxis.Visible = 'off';
        end
        
        % Adjust fontsize
        set(gca,'FontSize',20);
        
        % Set bar colours
        for b = 1:length(bar_obj)
            bar_obj(b).FaceColor = bar_colors(b,:);
        end
        
        % Add legend onto first image (second in experiment list)
        if k == 2
            leg_obj = legend({'Experiments','Param-free','Calibrated'},'FontSize',24);
            leg_obj.Position = [0.17, 0.84, 0.05, 0.035];
            leg_obj.Box = 'off';
        end
        
        % Label discrepancies if requested
        if label_discrepancies
            
            % Loop over each set of predictions
            for p = 1:length(predictions)
                
                % Add markers to the graph for discrepancies
                xtips = bar_obj(p+1).XEndPoints;
                ytips = bar_obj(p+1).YEndPoints;
                labels = cell(1,length(ind));
                labels(discrep_flags{p}(ind)) = {'x'};
                text(xtips,ytips,labels,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',24,'color',bar_colors(p+1,:))
                
            end
            
        end
        
        % Turn off the box, ticks, and thicken axis
        set(gca,'box','off','LineWidth',2,'TickLength',[0 0]);
        
        % Insert a dotted line if this plot contains beginning of observations
        % of the second component
        if ismember(twocomp_start,ind)
            
            % Calculate where the line would go on this plot (i.e. remove all
            % whole parts of plot_max) and put it between bars with a +0.5
            divide_pos = twocomp_start - floor(twocomp_start/plot_max)*plot_max + 0.5;
            line([divide_pos divide_pos],[-1 1.5*max(exp_y)], 'LineWidth',2, 'LineStyle','--','Color',[0 0 0]);
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y_trans = robustLinearTransform(y, z)
% This function returns a y that has been transformed from the input y to
% try to get closer to the value of z. This function assumes a simple
% linear regression of form
%      y_shift = beta_0 + beta_1 y
% but uses a more robust method to select beta_0 and beta_1, the Theil-Sen
% estimator. This method selects the slope using the median of all slopes
% between points and then selects an optimal intercept based on this slope

% Determine number of points provided
Np = length(y);

% Calculate slopes between all pairs of points (y,z)
m = zeros(Np*(Np-1)/2,1);
c = 0;
for i = 1:Np-1
    for j = i+1:Np
        
        % Calculate and store slope for this pair of points
        c = c + 1;
        m(c) = ( z(j) - z(i) ) / ( y(j) - y(i) );
        
    end
end

% Remove all NaN's from the slope calculation
m = m(~isnan(m));

% Find the median slope
beta_1 = median(m);

% Check to ensure slope is positive (rank-preserving)
if beta_1 < 0
    warning('A negative slope was returned, a different transformation (or none) will need to be applied.');
end

% Find the median residual to set the intercept
beta_0 = median( z - beta_1 * y );

% Use this transform to get new data
y_trans = beta_0 + beta_1 * y;