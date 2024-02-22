function result = runAllExperiments(theta, experiments)
% This function runs all the experiments for a given set of parameters and
% stores the results in a single cell array

% Use the model to predict results of all experiments
result = cell(1,length(experiments));
for k = 1:length(experiments) 
    result{k} = runExperiments(theta,@DunSS_paramODE,experiments{k});
end
