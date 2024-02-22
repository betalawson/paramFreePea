function RUNNER_DunABC


%%% Initial specifications

% Specify prior - uniform in the logs of the parameters
prior = priorUniform( [-3*ones(1,6), 1e-3, 1e-3], [3*ones(1,6),1,1] );


%%% Load in the data and prepare functions

% Read the information about experiments and their observations from data
[experiments, observations] = DunExperiments();

% Define the model function
f_model = @(theta) runAllExperiments([10.^theta(1:6),theta(7:8)], experiments);
% Use all model outputs as a summary
f_summaries = @(y) y;
% Define the discrepancy function
f_discrepancy = @(y, data) experimentDiscrepancy(y, data);


%%% Run the SMC

% Define options
options = struct('jumpType', 'MVN', 'Nparts', 5000, 'discrepancy', f_discrepancy, 'verbose', true);
% Call performSMC function
[particles, options, diagnostics] = performSMCABC( f_model, f_summaries, observations, prior, options );
% Save the results
save(['test_output','.mat'], 'particles', 'diagnostics', 'f_model', 'f_summaries', 'observations', 'prior', 'options');

end