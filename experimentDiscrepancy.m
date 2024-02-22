function [D, Ds] = experimentDiscrepancy(y,data)
% This function compares the model output in the input array 'y' and the
% experimental observations in the 'data' array. An optional second output
% provides the functionality of viewing individual discrepancies from each
% set of experimental data

% Initialise discrepancy
Ds = zeros(1,length(data));

% Specify whether to use cuttings or Kendall tau distance. The latter does
% not penalise ties
use_cuttings = true;

% Add discrepancy
for k = 1:length(data)
    
    % Gather together current data and eliminate non-observations
    d_here = data{k}(:);
    d_here(isnan(d_here)) = [];
    
    % Do the same for the model output
    y_here = y{k}(:);
    y_here(isnan(y_here)) = [];
    
    % If the two lists do not match in length, assign infinite discrepancy
    % (this happens when a steady state is not reached)
    if length(d_here) ~= length(y_here)
        Ds(k) = Inf;
    % Otherwise, discrepancy is total mismatch in ranking between the two
    else
        if use_cuttings
            Ds(k) = sum(findMismatches(y_here, d_here)) / length(y_here);
        else
            Ds(k) = kendallDist(d_here, y_here);
        end
    end

end

% Overall discrepancy is the mean discrepancy across experiments
D = mean(Ds);