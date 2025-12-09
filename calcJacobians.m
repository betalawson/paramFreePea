function [Jvar, Jparam] = calcJacobians(X, params, phi_s, phi_r)
% This function calculates the Jacobians (with respect to the model
% parameters and with respect to the state variables) of the parameterised
% version of the Dun et al. model for pea branching. The user inputs a 
% vector containing the values of the system's parameters together with the
% switching variables that specify the plant's genotype in its different
% compartments

% Read out the numbers of components from the genotypic switch inputs
Ns = size(phi_s,2);
Nr = size(phi_r,2);

% Count the number of variables - not including cytokinin
N_var = length(X) - Nr;

% Prepare the indexing used to refer to the different components of
% variable vector X
loc = 0;
FS_loc = loc + (1:Ns);
loc = loc + Ns;
fs_loc = loc + (1:Nr);
loc = loc + Nr;
SL_loc = loc + (1:Ns);
loc = loc + Ns;
sl_loc = loc + (1:Nr);
loc = loc + Nr;
I_loc = loc + (1:Ns);

% Read out these variables from their locations
FS = X(FS_loc);
fs = X(fs_loc);
SL = X(SL_loc);
sl = X(sl_loc);
I = X(I_loc);

% Read out the individual genotype variables
phi_FS = phi_s(1,:);
phi_SL = phi_s(2,:);
phi_I = phi_s(3,:);
phi_fs = phi_r(1,:);
phi_sl = phi_r(2,:);


%%%%%%%%%%%%%%%%%%%%%%%%% JACOBIAN FOR VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the variable Jacobian - 3 variables for shoot compartments, 
% two variables for each root compartment
Jvar = zeros(N_var,N_var);

% Derivatives for the equation for FS
Jvar(FS_loc, FS_loc) = diag(-1*ones(Ns,1));
Jvar(FS_loc, I_loc) = diag(params(1) .* 2 .* -params(2)./(params(2)+I).^2);
% Derivatives for the equation for fs
Jvar(fs_loc, fs_loc) = diag(-1*ones(Nr,1));
Jvar(fs_loc, FS_loc) = params(7) / Ns;
% Derivatives for the equation for SL
Jvar(SL_loc, SL_loc) = diag(-1*ones(Ns,1));
Jvar(SL_loc, FS_loc) = diag(2 * params(4) * FS .* phi_SL);
Jvar(SL_loc, sl_loc) = params(8) / Nr;
% Derivatives for the equation for sl
Jvar(sl_loc, sl_loc) = diag(-1*ones(Nr,1));
Jvar(sl_loc, fs_loc) = diag(2 * params(5) * fs .* phi_SL);
% Derivatives for the equation for I
Jvar(I_loc, I_loc) = diag(-1*ones(Ns,1));
Jvar(I_loc, SL_loc) = params(6) * phi_I;


%%%%%%%%%%%%%%%%%%%%%%%%% JACOBIAN FOR PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the parameter Jacobian - 3 variables for shoot compartments, 
% two variables for each root compartment, eight parameters
Jparam = zeros(N_var,8);

% Initialise tracking location
loc = 0;
% Derivatives with respect to parameters of equations for FS
Jparam(loc+(1:Ns),1) = 2 * params(2) ./ (params(2) + I) .* phi_FS;
Jparam(loc+(1:Ns),2) = 2 * params(1) .* I ./ (params(2) + I).^2 .* phi_FS;
loc = loc + Ns;
% Derivatives with respect to parameters of equations for fs
Jparam(loc+(1:Nr),3) = phi_fs;
Jparam(loc+(1:Nr),7) = 1/Ns * sum(FS);
loc = loc + Nr;
% Derivatives with respect to parameters of equations for SL
Jparam(loc+(1:Ns),4) = FS.^2 .* phi_SL;
Jparam(loc+(1:Ns),8) = 1/Nr * sum(sl);
loc = loc + Ns;
% Derivatives with respect to parameters of equations for sl
Jparam(loc+(1:Nr),5) = fs.^2 .* phi_sl;
loc = loc + Nr;
% Derivatives with respect to parameters of equations for I
Jparam(loc+(1:Ns),6) = SL .* phi_I;



end

