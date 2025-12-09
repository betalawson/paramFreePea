function FIGURE6_posteriorScatters()
% This function creates the figure used to show how the lower-discrepancy
% models are arranged throughout the parameter space.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the particles from the final output
load('final_particles.mat','particles','prior');

% Define the names of the variables
theta_names = {'log_{10} \lambda_{FS}','log_{10} K_{I}','log_{10} \lambda_{fs}','log_{10} \lambda_{SL}','log_{10} \lambda_{sl}','log_{10} \lambda_{I}','\alpha_{down}','\alpha_{up}'};
% Define the colours that mark the start and end of the Discrepancy scale
% NOTE: a linear interpolation in RGB space between these points is used
lowD_clr = [0.0 0.0 0.05];
highD_clr = [0.9 0.9 0.95];
% Specify how many rows and columns to use in plots
Nrows = 3;
Ncols = 5;
% Specify whether to plot particles in order of reducing discrepancy on
% scatterplots
sort_discrep = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% READ OUT DATA

% Raw particle data
thetas = getProperty(particles, 'theta');
Ds = getProperty(particles,'D');
Nthetas = size(thetas,2);
minD = min(Ds);
maxD = max(Ds);

% Sort particles by discrepancy so that best particles appear on top in
% scatterplots
if sort_discrep
    [Ds, I] = sort(Ds,'descend');
    thetas = thetas(I,:);
end


%%% PARAMETER BIVARIATE SCATTERS

% Figure setup
N_subfigs_per_fig = Ncols * Nrows;
Nsubfigs = Nthetas * (Nthetas - 1) / 2;
Nfigs = ceil(Nsubfigs / N_subfigs_per_fig);

% Create a colormap using the two colours specified
clrmap = [ linspace(lowD_clr(1), highD_clr(1), 50)', linspace(lowD_clr(2), highD_clr(2), 50)', linspace(lowD_clr(3), highD_clr(3), 50)' ];

% Create figure for parameter bivariate plots if it doesn't exist
for k = 1:Nfigs
    figObj = findobj('type','Figure','Name',['Bivariate Scatters ',num2str(k)]);
    if isempty(figObj)
        figs{k} = figure('units','Normalized','OuterPosition',[0 0 0.863 1], 'Name', ['Bivariate Scatters ',num2str(k)]);
    else
        figs{k} = figObj;
    end
end

% Initialise loop variables
c = 0;
cur_fig = 1;
% Move to current figure and clear it, then prepare axes
figure(figs{1});
clf;
tiles_obj = tiledlayout(Nrows,Ncols);
tiles_obj.Padding = 'compact';
tiles_obj.TileSpacing = 'compact';

% Loop over parameter combinations
for i = 1:Nthetas-1
    for j = i+1:Nthetas

        % Increment counter
        c = c + 1;
        
        % Move to new figure if counter exceeds limit
        if c > N_subfigs_per_fig
            c = 1;
            cur_fig = cur_fig + 1;
            figure(figs{cur_fig});
            clf;
            tiles_obj = tiledlayout(Nrows,Ncols);
            tiles_obj.Padding = 'compact';
            tiles_obj.TileSpacing = 'compact';
            nexttile;
            hold on;
        else
            nexttile;
            hold on;
        end
        
        % Plot scatter, with colour indicating discrepancy
        scatter(thetas(:,i), thetas(:,j), 15, Ds, 'filled');
        % Highlight best models in a separate color too
        scatter(thetas(Ds==minD,i), thetas(Ds==minD,j), 45, [1.0 0.7 0.3], 'filled');
        colormap(clrmap);
        xlabel(theta_names{i}, 'FontSize', 16);
        ylabel(theta_names{j}, 'FontSize', 16);
        set(gca, 'FontSize', 16);
        caxis([minD, maxD]);
        axis square;
        xlim([prior.theta_min(i), prior.theta_max(i)]);
        ylim([prior.theta_min(j), prior.theta_max(j)]);
        
    end
end

% Append a colorbar to the whole figure
clbar_obj = colorbar;
clbar_obj.Position = [0.62 0.175 0.33 0.035];
clbar_obj.Orientation = 'horizontal';
% Add text label above this
annotation('textbox',[0.63, 0.20, 0.13, 0.05], 'String','Discrepancy, \rho','EdgeColor','none','FontSize',18);


end