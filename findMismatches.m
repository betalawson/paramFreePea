function [discrep_flags, categories] = findMismatches(model_data, exp_data)
% This function looks through the model data, and determines the optimal
% way to 'pool' the model predictions, while preserving their rank
% ordering, such as to minimise discrepancies between the pooled
% experimental results and the model results

% Count number of experiments in this data
N = length(exp_data);

% Initialise all points as not a discrepancy to begin with
discrep_flags = false(1,N);

% Loop over all different values in the experimental data, and find best
% cuts in the model data values that can be made for each (this assumes the
% experimental data is already discretised into 1,2,3,...)
for m = 1:max(exp_data)-1
              
    % Find maximum model value for this experiment and minimum for next
    model_max = max(model_data(exp_data <= m));
    model_min = min(model_data(exp_data > m));
                
    % Determine discrepancies defined by two best possible cutoffs
    discreps1 = false(1,N);
    discreps2 = false(1,N);
    % Discrepancies according to cutoff 1
    if ~isempty(model_max)
        discreps1(exp_data > m & model_data < model_max + 1e-10) = true;
    % Or if this level doesn't exist in experiment, ignore it (by maximising discrepancy)
    else
        discreps1 = true(1,N);
    end
    % Discrepancies according to cutoff 2
    if ~isempty(model_min)
        discreps2(exp_data <= m & model_data > model_min - 1e-10) = true;
    % Or if this level doesn't exist in experiment, ignore it (by maximising discrepancy)
    else
        discreps2 = true(1,N);
    end
    
    % Use whichever cutoff is superior to determine discrepancies
    if sum(discreps1) < sum(discreps2)
        cutoffs(m) = model_max + 1e-10;
        discrep_flags(discreps1) = true;
    else
        cutoffs(m) = model_min - 1e-10;
        discrep_flags(discreps2) = true;
    end
    
end

% Apply the cutoffs found to classify the model data - work in descending
% order so that less than conditions can be simply applied in succession
categories = max(exp_data) * ones(size(model_data));
for m = max(exp_data)-1:-1:1
    categories( model_data < cutoffs(m) ) = m;
end