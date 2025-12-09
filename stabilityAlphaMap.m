function stabilityAlphaMap(varargin)
% This function takes as input either a genotype for a two-shoot, two-root
% plant, specified as a vector of 0's and 1's, or a genotype for any type 
% of plant, specified as two vectors phi_s and phi_r. Stability for a range
% of different alpha_up and alpha_down values is evaluated, and plotted

% Specify minimum alpha value
alpha_min = 1e-3;
% Specify number of points to evaluate in each direction
Npts = 101;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input handling
if nargin == 1
    
    G = varargin{1};
    phi_s = [G(1:2); G(3:4); G(5:6)];
    phi_r = [G(7:8); G(9:10)];
    
elseif nargin == 2
    
    phi_s = varargin{1};
    phi_r = varargin{2};
    
end

% Create a set of alpha values to run
alpha_up_vals = linspace(alpha_min, 1, Npts);
alpha_down_vals = linspace(alpha_min, 1, Npts);

% Loop over each and find its associated stability value
steady = zeros(Npts);
for i = 1:Npts
    
    alpha_up = alpha_up_vals(i);
    
    parfor j = 1:Npts
        steady(j,i) = DunStability([alpha_down_vals(j), alpha_up], phi_s, phi_r);
    end
    
end

%%% Visualise

% Prepare figure
figure('units','Normalized','OuterPosition',[0 0 1 1]);
% Plot results by padding plot with NaNs
[Aup, Adown] = meshgrid(alpha_up_vals, alpha_down_vals);
pcolor(Aup, Adown, steady);

% Clean up plot
axis square;
xlabel('\alpha_{up}','FontSize',24);
ylabel('\alpha_{down}','FontSize',24);
caxis([0 1]);
colormap([[1 0 0]; [1 1 1]]);
shading flat;
set(gca,'FontSize',24);

