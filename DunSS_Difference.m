function [SS, converged] = DunSS_Difference(alphas, phi_s, phi_r, delayed, provided_options)
% This function iterates the Dun et al. original difference equations to
% steady state (or a maximum number of iterations), and outputs both the
% equilibrium state reached (or NaNs if one is not reached) together with a
% flag that indicates whether a steady state was apparently reached. The
% user can supply a flag 'delayed' that specifies which version of the
% difference equation model to use.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify default options (can be overwritten by providing optional input
% argument 'options'
default_options.max_iters = 1e5;
default_options.SS_abstol = 1e-3;
default_options.SS_period = 5;
default_options.visualise = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume we are using the non-delayed version of the model if not specified
if nargin < 4
    delayed = false;
end

% Apply any user-provided options
options = default_options;
if nargin >= 5
    option_list = fieldnames(default_options);
    for k = 1:length(option_list)
        herename = option_list{k};
        if isfield(provided_options,herename)
            options.(herename) = provided_options.(herename);
        end
    end
end

% Read out how many compartments in shoot/root there are
Ns = size(phi_s,2);
Nr = size(phi_r,2);

% Initialise at a starting point shifted off of unity - dimension depends
% on whether or not old values also need to be stored
if delayed
    X_cur = 1.15 * ones( 1, (Ns+Nr)*4);
else
    X_cur = 1.15 * ones( 1, (Ns+Nr)*3);
end

% Initialise iteration count
iters = 0;
SS_iters = 0;

% Initialise at a state of not finding a steady state
converged = false;

% Prepare vector of all values for visualisation if doing so
if options.visualise
    Xv = NaN(options.max_iters+1,length(X_cur));
    Xv(1,:) = X_cur;
end

% Loop until a steady state is found, or max iterations reached
looping = true;
while looping
    
    % Increase iteration count
    iters = iters + 1;
    % Update old state
    X_old = X_cur;
    
    % Advance the difference equation
    X_cur = DunDifference(X_cur,alphas,phi_s,phi_r,delayed);
    
    % Check current state against previous
    if all( abs(X_cur - X_old) < options.SS_abstol )
        SS_iters = SS_iters + 1;
    else
        SS_iters = 0;
    end
    
    % Check if steady state has been reached
    if SS_iters >= options.SS_period
        looping = false;
        converged = true;
    end
    % Check for maximum number of iterations
    if iters >= options.max_iters
        looping = false;
    end
    
    % If visualising, store this next concentration in vector
    if options.visualise
        Xv(iters+1,:) = X_cur;
    end
    
end

% Initialise output to NaN's of appropriate size
SS = NaN(5,max([Ns,Nr]));
% Overwrite this with the real result if a steady state was reached
if converged
    loc = 0;
    SS(1,:) = X_cur(loc+(1:Ns));
    loc = loc + Ns;
    SS(2,:) = X_cur(loc+(1:Nr));
    loc = loc + Nr;
    SS(3,:) = X_cur(loc+(1:Ns));
    loc = loc + Ns;
    SS(4,:) = X_cur(loc+(1:Nr));
    loc = loc + Nr;
    SS(5,:) = X_cur(loc+(1:Ns));
    loc = loc + Ns;
    SS(6,:) = X_cur(loc+(1:Nr));
end

% Plot results if visualising
if options.visualise
    plot(Xv);
end

end