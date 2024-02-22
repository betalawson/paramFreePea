function options = DEFINE_SMCdefaults
% This function is where the user can set the options for the SMC

options.Nparts = 5000;                  % Number of particles
options.alpha = 0.75;                   % Keep ratio
options.minMCMCsteps = 10;              % Minimum number of MCMC steps (number to perform before estimating number of steps to take)
options.maxMCMCsteps = 1000;            % Maximum number of MCMC steps

options.jumpType = 'MVN';               % Type of jumping distribution to use
options.meanFlip = 2;                   % Maximum number of bits to flip for binary updates

%%%%% SMC-ABC specific %%%%%
options.discrepancy = 'mahalanobis';    % Discrepancy function - may provide a function @(S) that returns the discrepancy for a set of summary statistics, or may specify 'euclid', 'weighted', 'mahalanobis', 'SSnonans'
options.mutateAll = 'true';             % Specifies whether all particles will have their positions updated by the mutation step (used in SMC-ABC)
options.terminate_D = false;            % Specifies whether to terminate SMC-ABC when the target discrepancy (set below) is reached
options.terminate_acceptance = true;    % Specifies whether to terminate SMC-ABC when the maximum number of MCMC steps is reached
options.terminate_progress = true;      % Specifies whether to terminate SMC-ABC when progress (judged by evolving threshold of discrepany thresholds) has stalled
options.D_target = -Inf;                % Target discrepancy that triggers exit from the SMC-ABC routine if 'terminate_D' flag set
options.stall_tol = 1e-4;               % The minimum level of decrease required for an iteration to count as not stalled
options.stall_iterations = 5;                % The number of non-improving iterations that triggers termination

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.verbose = true;                     % Display output while running
options.visualise = true;                   % Display graphics while running
options.visFunc = @defaultVisualisation;    % Function to use for graphical display

end