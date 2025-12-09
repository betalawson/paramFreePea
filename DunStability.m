function stable = DunStability(alphas, phi_s, phi_r)
% This function considers the viability and the stability of the fixed
% points of the first (simpler) Dun model
%
% The inputs are the proportion of useful transport upwards and downwards
% (alpha_up and alpha_down), and the state of the genetic switches, g

% Specify the tolerance for equilibria (root mean square)
tol = 5e-4;
% Specify initial number of iterations to try in finding a steady state
init_iters = 500;
% Specify maximum number of fminsearch iterations
max_minsearch_iters = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First try the iterator with a small number of iterations to check for
% stability that is reached quickly
SS = DunSS_Difference(alphas, phi_s, phi_r, true, struct('max_iters',init_iters));
if ~all(isnan(SS),'all')
    stable = true;
    stability_found = true;
else
    stability_found = false;
end

if ~stability_found
    
    % Read out how many compartments in shoot/root there are
    Ns = size(phi_s,2);
    Nr = size(phi_r,2);
    
    % Calculate number of total variables
    Nx = 4 * (Ns+Nr);
    
    % Create a function that evalutes the difference equation RHS
    DiffRHS = @(X) DunDifference(X,alphas,phi_s,phi_r,true);
    
    % Find the steady state point
    options = optimset('MaxFunEvals',1e4,'Display','off');
    Xstar = fminsearch(@(X) rms(DiffRHS(X) - X),1.15*ones(1,Nx),options);
    % Repeatedly apply fminsearch re-initialised as this works better
    err = rms(DiffRHS(Xstar) - Xstar);
    minsearch_iters = 1;
    while err > tol && minsearch_iters < max_minsearch_iters
        minsearch_iters = minsearch_iters + 1;
        Xstar = fminsearch(@(X) rms(DiffRHS(X) - X),Xstar,options);
        err = rms(DiffRHS(Xstar) - Xstar);
    end
    
    % Check the minimum found to see if it really solves the equations
    err = rms(DiffRHS(Xstar) - Xstar);
    if err < tol
        
        % Numerically estimate the Jacobian at this point
        J = zeros(Nx);
        for i = 1:Nx
            
            % Create a perturbation in one direction
            pert = [zeros(1,i-1),sqrt(eps),zeros(1,Nx-i)];
            % Numerically evaluate derivative here
            J(:,i) = (DiffRHS(Xstar + pert) - DiffRHS(Xstar - pert)) / (2 * sqrt(eps));
            
        end
        
        % Use Jacobian eigenvalues to evaluate stability, - if they are
        % conclusive
        max_modulus = max(abs(eig(J)));
        if max_modulus < 0.99
            stable = true;
            stability_found = true;
        elseif max_modulus > 1.01
            stable = false;
            stability_found = true;
        end
        
    end
    
end

if ~stability_found
        
    % Run the steady state checker code instead
    SS = DunSS_Difference(alphas, phi_s, phi_r, true);
    if ~all(isnan(SS),'all')
        stable = true;
    else
        stable = false;
    end
        
end
    
end