function FIGURE_SMCABC_discrepHistory(diagnostics)
% This function uses the provided set of diagnostics from an SMC-ABC run to
% visualise the history of particle discrepancies over "time" (iterations
% of the sampling method)

% Calculate the base model discrepancy
[experiments, observations] = DunExperiments();
base_params = [1, 1, 1, 1, 1, 1, 0.5, 0.5];
base_output = runAllExperiments(base_params, experiments);
base_D = experimentDiscrepancy(base_output, observations);

% Extract the particle discrepancy history from the provided structure
particles_history = diagnostics.particles_history;
N_history = length(particles_history);
for k = 1:N_history
    % Grab out particles for this iteration
    particles_here = particles_history{k};
    % Extract discrepancy values for this particle set
    part_Ds_history(:,k) = getProperty(particles_here,'D');
end

% Calculate the range of values
minDs = min(part_Ds_history);
maxDs = max(part_Ds_history);
% Also extract the median, and 95% credible interval
minQ = quantile(part_Ds_history,0.025);
maxQ = quantile(part_Ds_history,0.975);
medQ = median(part_Ds_history);

% Prepare extra data used for filling in the curve
% This has some extra dummy elements added just to remove the patch
% boundaries on left and right (hacky, yes)
x = (1:N_history) - 1;
xflip = [[-1,x,N_history+1], fliplr([-1,x,N_history+1])];
inBetween = [[minQ(1),minQ,minQ(end)], fliplr([maxQ(1),maxQ,maxQ(end)])];

% Now, plot all of these quantities together
figure; hold on;

% Plot the 95% credible interval region
patch_obj = fill(xflip, inBetween, [0.8 0.8 1.0]);
patch_obj.EdgeColor = [0.6 0.6 0.8];
patch_obj.LineWidth = 2.5;
% Plot the minimum and maximum values using softer-coloured dotted lines
plot(x,minDs,'--.','LineWidth',2.5,'MarkerSize',30,'Color',[0.6 0.6 0.8]);
plot(x,maxDs,'--.','LineWidth',2.5,'MarkerSize',30,'Color',[0.6 0.6 0.8]);
% Plot the median value
plot(x,medQ,'.-','LineWidth',4,'MarkerSize',40,'Color',[0.0, 0.0, 0.9],'MarkerEdgeColor',[0.0, 0.0, 0.9]);
% Plot the base discrepancy value
plot(x,base_D*ones(1,N_history),'k--','LineWidth',4);

% Add details to plot and fix axes
xlim([0 N_history-1]);
xlabel('SMC-ABC Iterations');
ylabel('Discrepancy with Qualitative Data');
set(gca,'FontSize',24);


