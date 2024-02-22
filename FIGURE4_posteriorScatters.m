function FIGURE4_posteriorScatters()
% This function creates the figure used to show how the lower-discrepancy
% models are arranged throughout the parameter space.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the particles from the final output
load('final_particles.mat','particles');

% Define the names of the variables
theta_names = {'log_{10} \lambda_{FS}','log_{10} K_{I}','log_{10} \lambda_{fs}','log_{10} \lambda_{SL}','log_{10} \lambda_{sl}','log_{10} \lambda_{I}','\alpha_{down}','\alpha_{up}'};
% Define the colours that mark the start and end of the Discrepancy scale
% NOTE: a linear interpolation in RGB space between these points is used
lowD_clr = [0.0 0.0 0.05];
highD_clr = [0.9 0.9 0.95];
% Specify how many rows and columns to use in plots
Nrows = 4;
Ncols = 8;
% Specify whether to plot particles in order of reducing discrepancy on
% scatterplots
sort_discrep = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% READ OUT DATA

% Raw particle data
thetas = getProperty(particles, 'theta');
Ds = getProperty(particles,'D');
Nthetas = size(thetas,2);

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
        figs{k} = figure('units','Normalized','OuterPosition',[0 0 1 1], 'Name', ['Bivariate Scatters ',num2str(k)]);
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
ax = createAxes( Nsubfigs, Nrows, Ncols, [], true );

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
            ax = createAxes( Nsubfigs, Nrows, Ncols, [], true );
        end
        
        % Plot scatter, with colour indicating discrepancy
        scatter(ax{c}, thetas(:,i), thetas(:,j), 20, Ds, 'filled');
        colormap(clrmap);
        xlabel(ax{c}, theta_names{i}, 'FontSize', 16);
        ylabel(ax{c}, theta_names{j}, 'FontSize', 16);
        set(ax{c}, 'FontSize', 16);
        
    end
end

% Append a colorbar to the whole figure
clbar_obj = colorbar;
clbar_obj.Position = [0.552 0.19 0.35 0.035];
clbar_obj.Orientation = 'horizontal';
% Add text label above this
annotation('textbox',[0.577, 0.22, 0.15, 0.05], 'String','Discrepancy, \rho','EdgeColor','none','FontSize',18);


end