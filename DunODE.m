function dX = DunODE(t, X, alphas, phi_s, phi_r)
% This function evaluates the right hand side of the ODE model associated
% with the Dun et al. difference equation, using the transfer parameters
% provided in vector 'alphas' and genotypic switch settings taking
% values {0,1}:
%  phi_s: 3xNs, Ns the number of shoot components
%         (each column specifies [phi_FS; phi_SL; phi_I] for one component)
%  phi_r: 2xNr, Nr the number of root components
%         (each column specifies [phi_fs; phi_sl] for one component)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise the output vector
dX = zeros(size(X));

% Read out the numbers of components from the genotypic switch inputs
Ns = size(phi_s,2);
Nr = size(phi_r,2);

% Read out the genotypic switch settings and make them column vectors
phi_FS = phi_s(1,:)';
phi_SL = phi_s(2,:)';
phi_I = phi_s(3,:)';
phi_fs = phi_r(1,:)';
phi_sl = phi_r(2,:)';

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
ck = X(loc+(1:Nr));

% Use the Dun et al. model equations to derive expressions for dX/dt
dFS = 2  ./ (1 + I) .* phi_FS - FS;
dfs = phi_fs + alphas(1) / Ns * sum(FS) - fs;
dSL = FS.^2 .* phi_SL + alphas(2) / Nr * sum(sl) - SL;
dsl = fs.^2 .* phi_sl - sl;
dI = SL .* phi_I - I;
dck = 2 ./ (1 + fs) - ck;

% Re-store the individual components of dX in dX
loc = 0;
dX(loc+(1:Ns)) = dFS;
loc = loc + Ns;
dX(loc+(1:Nr)) = dfs;
loc = loc + Nr;
dX(loc+(1:Ns)) = dSL;
loc = loc + Ns;
dX(loc+(1:Nr)) = dsl;
loc = loc + Nr;
dX(loc+(1:Ns)) = dI;
loc = loc + Ns;
dX(loc+(1:Nr)) = dck;

end