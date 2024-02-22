function stabilityManyPlants()

% Specify how many plants (values for alphas) to run per genotype
N_plants_per_genotype = 1000;
% Specify minimum value of alpha to consider
alpha_min = 1e-3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% TWO SHOOT, TWO ROOT PLANTS

% Create all combinations of genotypic state for two-shoot, one-root plants
[X1, X2, X3, X4, X5, X6, X7, X8, X9, X10] = ndgrid([0,1]);
% Gather this all into one large genotype matrix
G = [X1(:), X2(:), X3(:), X4(:), X5(:), X6(:), X7(:), X8(:), X9(:), X10(:)];
% Specify the paired variables (corresponding to matched quantities in
% multiple components). First set are shoot, second set are root.
pairs = {[1, 2; 3, 4; 5, 6], [7, 8; 9, 10]};
% Cut out any redundant genotypes (matched when swapping across all pairs
% in a set)
G = removeRedundantCombinations(G,pairs);

phiFS1_states = G(:,1);
phiFS2_states = G(:,2);
phiSL1_states = G(:,3);
phiSL2_states = G(:,4);
phiI1_states = G(:,5);
phiI2_states = G(:,6);
phifs1_states = G(:,7);
phifs2_states = G(:,8);
phisl1_states = G(:,9);
phisl2_states = G(:,10);

% Count number of genotypes for this graft mix
N_genotypes = length(phiFS1_states);

% Loop over each, determining stability of the plant with that genotype
% across many different values for the two parameters (alpha's)
s2r2_steady = zeros(1,N_genotypes);
for k = 1:N_genotypes
    
    fprintf('Working on genotype %g out of %g combinations...',k,N_genotypes);
    
    % Compile genotype information into shoot and root phenotype info
    phi_s = [phiFS1_states(k), phiFS2_states(k);
        phiSL1_states(k), phiSL2_states(k);
        phiI1_states(k) , phiI2_states(k)   ];
    phi_r = [phifs1_states(k), phifs2_states(k);
        phisl1_states(k), phisl2_states(k)  ];
    
    % If this is a known special case, process it separately (helps speed)
    if all(phi_s(:,1)) && ~any(phi_r(2,:))
        s2r2_steady(k) = 0;    % Eigenvalue of exactly one
    elseif ~phi_s(1,1) || ~phi_s(3,1) || (~phi_s(2,1) && ~phi_r(2,1))
        s2r2_steady(k) = 1;    % Non-functional plant
    else
        
        % Select randomly many combinations of (alpha1, alpha2) to trial
        alphas = alpha_min + (1-alpha_min)*rand(N_plants_per_genotype-4,2);
        
        % Append the most extreme cases to ensure they are tested
        alphas = [alphas; alpha_min, alpha_min; alpha_min, 1; 1, alpha_min; 1, 1];
        
        % Repeatedly create random plants with this genotype and evaluate
        steady_flags = false(1,N_plants_per_genotype);
                
        parfor n = 1:N_plants_per_genotype
            
            % Grab out current alpha values
            alphas_here = alphas(n,:);
            % Check stability for these alpha values
            steady_flags(n) = DunStability(alphas_here, phi_s, phi_r);
            
        end
        
        % Store the proportion of steady plants for this genotype
        s2r2_steady(k) = mean(steady_flags);
    end
   
    fprintf(' proportion of stable alpha values: %g%%\n',s2r2_steady(k)*100);
    
end

s2r2_genotypes = [phiFS1_states, phiFS2_states, phiSL1_states, phiSL2_states, phiI1_states, phiI2_states, phifs1_states, phifs2_states, phisl1_states, phisl2_states];

% Save all results
save('many_plants_stability.mat','s2r2_steady','s2r2_genotypes');

end

function G = removeRedundantCombinations(G,pairs)

% First, clean data by removing any duplicated combinations. Duplicated
% combinations are after switching the labels of all paired components (in
% shoot or root), the result is an already-recorded entry

% Create a shoot pair-swapped version of the genotype list
Gs = G;
Gs(:,pairs{1}(:,1)) = G(:,pairs{1}(:,2));
Gs(:,pairs{1}(:,2)) = G(:,pairs{1}(:,1));
% Create a root pair-swapped version of the genotype list
Gr = G;
Gr(:,pairs{2}(:,1)) = G(:,pairs{2}(:,2));
Gr(:,pairs{2}(:,2)) = G(:,pairs{2}(:,1));
% Create a root+shoot pair-swapped version of the genotype list
Gsr = Gs;
Gsr(:,pairs{2}(:,1)) = Gs(:,pairs{2}(:,2));
Gsr(:,pairs{2}(:,2)) = Gs(:,pairs{2}(:,1));
% Gather together all lists to search over
search_lists = {Gs, Gr, Gsr};

% Loop through the data, removing entries match to within a swap
looping = true;
loc = 0;
while looping
    
    % Move along vector
    loc = loc + 1;
    % Grab out current genotype
    Gc = G(loc,:);
        
    % Search for any genotypes coinciding in any of the swapped lists
    for i = 1:length(search_lists)
        % Find any matches, avoiding current location
        [found, matchloc] = ismember(Gc,search_lists{i}([1:loc-1,loc+1:end],:),'rows');
        % Correct matching location to account for the deleted element
        if matchloc >= loc
            matchloc = matchloc + 1;
        end
        % If a match is found, delete it
        if found
            G(matchloc,:) = [];
            % Corret current location in vector to account for deletion
            if matchloc <= loc
                loc = loc - 1;
            end
            % Delete from search list
            for j = 1:length(search_lists)
                search_lists{j}(matchloc,:) = [];
            end
        end
    end
    % Terminate loop if end of vector is hit
    if loc == size(G,1)
        looping = false;
    end
    
end

end