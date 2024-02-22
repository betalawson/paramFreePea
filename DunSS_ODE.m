function [SS, converged] = DunSS_ODE(alphas, phi_s, phi_r)
% This function iterates the Dun ODE to steady state (or a maximum
% specified time), and outputs both the equilibrium state reached (or NaNs
% if one is not reached) together with a flag that indicates whether a
% steady state was apparently reached

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify the maximum end time
t_max = 1000;
t_segment = 5;

% Specify the tolerance used for steady state detection
SS_abstol = 1e-5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read out the number of shoot and root compartments
Ns = size(phi_s,2);
Nr = size(phi_r,2);

% Initialise at a starting point shifted off of unity
X_cur = 1.15 * ones( (Ns+Nr)*3, 1);

% Initialise remaining time at maximum time
t_remain = t_max;

% Initialise at a state of not finding a steady state
converged = false;

% Loop until a steady state is found, or max time value is reached
looping = true;
while looping
    
    % Simulate the ODE for the next segment of time
    [T,X] = ode15s( @(t,X) DunODE(t,X,alphas,phi_s,phi_r), [0 t_segment], X_cur);
    % Update time remaining and state of X
    t_remain = t_remain - t_segment;
    X_cur = X(end,:)';
    % Check for steady state across the whole segment
    dXdt = diff(X)./diff(T);
    if all( abs(dXdt) < SS_abstol, 'all')
        looping = false;
        converged = true;
    end
    % Check for time exceeded
    if t_remain <= 0
        looping = false;
    end

end

% Initialise output to NaN's of appropriate size
SS = NaN(6,max([Ns,Nr]));
% Overwrite this with the real result if a steady state was reached
if converged
    loc = 0;
    SS(1,1:Ns) = X_cur(loc+(1:Ns))';
    loc = loc + Ns;
    SS(2,1:Nr) = X_cur(loc+(1:Nr))';
    loc = loc + Nr;
    SS(3,1:Ns) = X_cur(loc+(1:Ns))';
    loc = loc + Ns;
    SS(4,1:Nr) = X_cur(loc+(1:Nr))';
    loc = loc + Nr;
    SS(5,1:Ns) = X_cur(loc+(1:Ns))';
    loc = loc + Ns;
    SS(6,1:Nr) = X_cur(loc+(1:Nr))';
end

end