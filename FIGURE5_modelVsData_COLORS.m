function FIGURE5_modelVsData_COLORS
% This function produces the figure presenting comparisons between the best
% particle's model, the baseline (parameter-free) model, and the data from
% experiments. Given the variability in scale of model predictions, these
% predictions are mapped onto values similar to the discrete values
% recorded for the experiments using a linear transformation that preserves
% rank and thus does not affect the discrepancy measure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify box colours for the plots. Wild type is white, upregulation is
% shades of blue and downregulation is shades of red
level_colors3 = [1.0, 0.4, 0.4;
                 1.0, 1.0, 1.0;
                 0.4, 0.4, 1.0];   
level_colors5 = [1.0, 0.4, 0.4;
                 1.0, 0.7, 0.7;
                 1.0, 1.0, 1.0;
                 0.7, 0.7, 1.0;
                 0.4, 0.4, 1.0];   

% List the species contained in the experiments
species = {'Cytokinin','{\it RMS1} (shoot)','{\it RMS1} (root)','Branching'};

% List the predictions being plotted
comparison_data = {'Experimental', 'Parameter-free', 'Parameterised'};

% Specify the number of boxes to plot on each figure (if not enough data
% present, "ghost boxes" will be plotted for sizing reasons)
N_plot = 43;

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
    
    % Write the experimental data as first row in a storage matrix
    results = exp_y(:)'; 
    % Loop over model predictions, classifying them and appending these
    for p = 1:length(predictions)
        
        % Read out the model's predictions for this experiment
        model_y = predictions{p}{k};
        
        % Apply same re-ordering applied to data to the model output
        model_y = model_y(I,:);
        
        % Grab out only the valid ones using the template
        model_y = model_y(template);
        
        % Classify model output and append classified results to matrix
        [~,categories] = findMismatches(model_y,exp_y);
        results(end+1,:) = categories(:)';
        
    end
    
    % Italicise all mutant names
    for r = 1:5
        exp_mutants = replace(exp_mutants,['rms',num2str(r)],['{\itrms',num2str(r),'}']);
        exp_mutants = replace(exp_mutants,['r',num2str(r)],['{\itr',num2str(r),'}']);
    end
    
    % For anything but branching inhibition, shift the results so that the
    % wild type is zero
    if ~strcmp(species{k},'Branching')
        results = results - (max(exp_y)+1)/2;
    else
        % For branching inhibition, move [1,3] data to [-2,0] then negate
        results = -(results - 3);
    end
        
    % Italicise all mutant names
    for r = 1:5
        exp_mutants = replace(exp_mutants,['rms',num2str(r)],['{\itrms',num2str(r),'}']);
        exp_mutants = replace(exp_mutants,['r',num2str(r)],['{\itr',num2str(r),'}']);
    end
    
    % Flip results so experiments are on top
    results = flipud(results);
    
    % Now loop over chunks of the data, ensuring all are plotted to the
    % same specified size above
    plotting_complete = false;
    loc = 0;
    while ~plotting_complete
    
        % Create a new figure
        fig_obj = figure('units','Normalized','OuterPosition',[0 0 1 1]);
        ax_obj = gca(fig_obj);
        
        % Check how much data there is to actually plot
        N_here = min([size(results,2)-loc, N_plot]);
        
        % Create results container with NaN data to be filled in
        results_here = nan(size(results,1),N_plot);
        results_here(:,1:N_here) = results(:,loc+(1:N_here));
        
        % Read out mutants for these results
        mutants_here = exp_mutants(loc+(1:N_here));
    
        % Plot using pcolor, after appending dummy data because pcolor
        results_here(end+1,:) = NaN;
        results_here(:,end+1) = NaN;
        pcolor(results_here);
        shading flat;
        box off;
    
        % Use colormap based on number of levels
        if max(results(:)) == 1
            colormap(level_colors3);
            caxis([-1, 1]);
        else
            colormap(level_colors5);
            caxis([-2, 2]);
        end
        
        % Make the plot squares instead of rectangles
        axis equal;
    
        % Trim graph to avoid ugly padding
        ylim(ax_obj, [1 length(predictions)+2]);
        
        % Provide ticks along x - each experiment's genotype
        xticks(ax_obj, (1:N_here) + 0.5 );
        xticklabels(ax_obj, mutants_here);
        ax_obj.XAxis.TickLabelRotation = 45;
    
        % Add y labels except for rms1-root, which sits beside rms1-shoot
        if ~strcmp(species{k},'rms1 (root)')
        
            % Provide ticks along y - the data being plotted on that row
            yticks(ax_obj, (1:length(predictions)+1) + 0.5 );
            % Flip list of plotted data to match flip applied to results matrix
            yticklabels(ax_obj, fliplr(comparison_data) );
            ax_obj.YAxis.FontSize = 20;
            
        else
            
            yticks(ax_obj, []);
            
        end
        
        % Move ticks to outside
        ax_obj.TickDir = 'out';
        % Thicken axis lines and remove Y-tick markers
        ax_obj.LineWidth = 2;
        ax_obj.YAxis.TickLength = [0 0];
                
        % Add dotted line and two titles if two components present
        if 0 < twocomp_start - loc && twocomp_start - loc < N_plot
           
            % Two component location is the value minus previously plotted
            twocomp_loc = twocomp_start - loc;
            
            % Add grey line
            hold(ax_obj,'on');
            line(ax_obj,[1 1] * twocomp_loc,[1-0.6 length(predictions)+2+0.6],'LineWidth',3,'Color',[0.4 0.4 0.4],'LineStyle','--'); 
            ax_obj.Clipping = 'off';
            % Add two titles
            text(ax_obj,twocomp_loc / 2 + 1, length(predictions)+2+0.6,species{k},'FontSize',24,'FontWeight','bold','HorizontalAlignment','center');
            text(ax_obj,twocomp_loc + (N_plot - twocomp_loc)/2 + 1, length(predictions)+2+0.6,'(Second Component)','FontSize',24,'FontWeight','bold','HorizontalAlignment','center');
            
        else
            
            % Add title and centre it above actual plotted data
            text(ax_obj,N_here/2 + 1, length(predictions)+2+0.6,species{k},'FontSize',24,'FontWeight','bold','HorizontalAlignment','center');
        
        end
                
        % Check loop, add to location if data remains
        if size(results,2) <= loc+N_plot
            plotting_complete = true;
        else
            loc = loc + N_plot;
        end
        
    end
        
end