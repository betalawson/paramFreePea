function X = DunDifference(X, alphas, phi_s, phi_r, delayed)
% This function evaluates the right hand side of the ODE model associated
% with the Dun et al. difference equation, using the transfer parameters
% provided in 2x1 vector 'alphas' and genotypic switch settings taking
% values {0,1}:
%  phi_s: 3xNs, Ns the number of shoot components
%         (each column specifies [phi_FS; phi_SL; phi_I] for one component)
%  phi_r: 2xNr, Nr the number of root components
%         (each column specifies [phi_fs; phi_sl] for one component)
% delayed: true/false specifying whether transport is delayed an iteration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assume no delay if the flag isn't provided
if nargin < 5
    delayed = false;
end

% Read out the numbers of components from the genotypic switch inputs
Ns = size(phi_s,2);
Nr = size(phi_r,2);

% Read out the genotypic switch settings
phi_FS = phi_s(1,:);
phi_SL = phi_s(2,:);
phi_I = phi_s(3,:);
phi_fs = phi_r(1,:);
phi_sl = phi_r(2,:);

% Read out from the dependent variable vector X the different equation
% components
loc = 0;
FS = X(loc+(1:Ns));
loc = loc + Ns;
fs = X(loc+(1:Nr));
loc = loc + Nr;
SL = X(loc+(1:Ns));
loc = loc + Ns;
sl = X(loc+(1:Nr));
loc = loc + Nr;
I = X(loc+(1:Ns));
loc = loc + Ns;
ck = X(loc+(1:Nr));  % cytokinin is purely a model output, hence warning

% Continue reading out old states if using the delayed equation
if delayed
    loc = loc + Nr;
    FS_old = X(loc+(1:Ns));
    loc = loc + Ns;
    sl_old = X(loc+(1:Nr));
end
    

%%% DUN ET AL. MODEL DEFINITION

% Shoot feedback signal
FS_new = 2 ./ (1 + I) .* phi_FS;

% Root feedback signal
if delayed
    fs_new = phi_fs + alphas(1) / Ns * sum(FS_old);
else
    fs_new = phi_fs + alphas(1) / Ns * sum(FS);
end

% Shoot strigolactone
if delayed
    SL_new = FS.^2 .* phi_SL + alphas(2) / Nr * sum(sl_old);
else
    SL_new = FS.^2 .* phi_SL + alphas(2) / Nr * sum(sl);
end

% Root strigolactone
sl_new = fs.^2 .* phi_sl;

% Shoot inhibitory factor
I_new = SL .* phi_I;

% Root cytokinin
ck_new = 2 ./ (1 + fs);


% Re-store the individual components of dX in dX
loc = 0;
X(loc+(1:Ns)) = FS_new;
loc = loc + Ns;
X(loc+(1:Nr)) = fs_new;
loc = loc + Nr;
X(loc+(1:Ns)) = SL_new;
loc = loc + Ns;
X(loc+(1:Nr)) = sl_new;
loc = loc + Nr;
X(loc+(1:Ns)) = I_new;
loc = loc + Ns;
X(loc+(1:Nr)) = ck_new;
% Also store delayed variables if using that model
if delayed
    loc = loc + Nr;
    X(loc+(1:Ns)) = FS;
    loc = loc + Ns;
    X(loc+(1:Nr)) = sl;
end

end