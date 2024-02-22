function manyPlantsSummariser()
% This function loads in the data for multiple plants and outputs it in a
% way that might help the user identify trends

% Load in the data
load('many_plants_stability.mat','s2r2_genotypes','s2r2_steady');

%%% TWO SHOOT TWO ROOT

% Prepare data
genodata = s2r2_genotypes;
steady_props = s2r2_steady;
genotype_names = {'FS1', 'FS2', 'SL1', 'SL2', 'I1', 'I2', 'fs1', 'fs2', 'sl1', 'sl2'};


% Text to guide user
fprintf('\n------------------------\n');
fprintf('TWO SHOOT TWO ROOT CASE\n');
fprintf('------------------------\n\n');

% Run code
genotypeSteadySummarise(genodata, steady_props, genotype_names, settings);

end

function genotypeSteadySummarise(G, p, names, settings)

% Count total number of unique genotypes
Ngen = size(G,1);

% Grab out stable, unstable, and ambiguous genotypes
Gs = G(p > 1 - 1e-10,:);
Gu = G(p < 1e-10,:);
Ga = G(p >= 1e-10 & p <= 1 - 1e-10,:);
pa = p(p >= 1e-10 & p <= 1 - 1e-10);

% First just summarise the numbers in each category
fprintf('Across all unique genotype combinations, %g%% of plants are stable, %g%% are unstable, and %g%% are alpha-dependent\n', size(Gs,1)/Ngen*100, size(Gu,1)/Ngen*100, size(Ga,1)/Ngen*100);
% Output proportions of state for each gene in the different categories
Gs_onprop = mean(Gs,1)
Gu_onprop = mean(Gu,1)
Ga_onprop = mean(Ga,1)

% Grab out plants that match in all components - already analysed in paper
single_comp = ( G(:,1) == G(:,2) ) & ( G(:,3) == G(:,4) ) & ( G(:,5) == G(:,6) ) & ( G(:,7) == G(:,8) ) & ( G(:,9) == G(:,10) );
N1 = sum(single_comp);

% Similarly analyse these
G1 = G(single_comp,:);
p1 = p(single_comp);
G1s = G1(p1 > 1 - 1e-10,:);
G1u = G1(p1 < 1e-10,:);
G1a = G1(p1 >= 1e-10 & p1 <= 1 - 1e-10,:);
fprintf('Across one-component genotype combinations, %g%% of plants are stable, %g%% are unstable, and %g%% are alpha-dependent\n', size(G1s,1)/N1*100, size(G1u,1)/N1*100, size(G1a,1)/N1*100);
% Output proportions of state for each gene in the different categories
G1s_onprop = mean(G1s,1)
G1u_onprop = mean(G1u,1)
G1a_onprop = mean(G1a,1)

% Remaining plants are the novel two-component plants
G2 = G(~single_comp,:);
p2 = p(~single_comp);
N2 = sum(~single_comp);
G2s = G2(p2 > 1 - 1e-10,:);
G2u = G2(p2 < 1e-10,:);
G2a = G2(p2 >= 1e-10 & p2 <= 1 - 1e-10,:);
fprintf('Across novel two-component genotype combinations, %g%% of plants are stable, %g%% are unstable, and %g%% are alpha-dependent\n', size(G2s,1)/N2*100, size(G2u,1)/N2*100, size(G2a,1)/N2*100);
% Output proportions of state for each gene in the different categories
G2s_onprop = mean(G2s,1)
G2u_onprop = mean(G2u,1)
G2a_onprop = mean(G2a,1)


keyboard

end




