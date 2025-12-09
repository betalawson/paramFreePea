function FIGURE3_stabilityMaps()
% This function re-creates the results shown in Figure 3, which displays
% the map of stability of steady state under some genotypes. Note that this
% is provided for reproducibility, several modifications were made to
% figures (adding textbox labels, fontsizes, axes widths, etc.) to improve
% graphical appearance that have not been included here

% Run the stability map generator for one example that shows
% alpha-dependent stability
stabilityAlphaMap([1, 1; 1, 0; 1, 0],[0, 0; 1, 0]);

% Run the specific example for which high values of alpha_up and alpha_down
% also cause the steady state to become unstable
stabilityAlphaMap([1, 1; 1, 0; 1, 1],[0; 1]);

end

