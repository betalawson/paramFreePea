function outputs = runExperiments(params, f_model, experiments)
% This function runs all the experiments implied by the provided struct and
% collects their outputs

% Read out how many experiments were provided
N_exp = length(experiments.phi_s);

% Storage initialisation
outputs = NaN(N_exp,2);

% Run each experiment in turn
for r = 1:length(experiments.phi_s)
    
    % Find steady state
    SS = f_model(params, experiments.phi_s{r}, experiments.phi_r{r});
    % Get result of this experiment
    result = experiments.f_measure{r}(SS);
    % Store this
    outputs(r,1:length(result)) = result;
    
end

end

