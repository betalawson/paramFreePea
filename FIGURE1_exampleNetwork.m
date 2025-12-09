function FIGURE1_exampleNetwork


%%% SIMULATION SPECIFICATIONS
bool_IC = [1,0,0,0];
ode_IC = [1,0,0.5,0];
N_gen = 16;
off_gen = 8;


%%% PLOT SPECIFICATIONS
plot_colors = [0.4, 0.4, 1.0;
               0.6, 0.0, 0.0;
               1.0, 0.6, 0.6;
               0.0, 0.0, 0.0];
           
species_names = {'MYC2', 'EPF2', 'EPFL9', 'ERECTA'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% BOOLEAN NETWORK MODEL

% Create the list of all possible states (binary numbers 0-15)
state_list = arrayfun( @(x) str2num(x), dec2bin(0:15));
% Implementation of transition rules - default to states remaining fixed
new_state_list = state_list;
% First rule - state of MYC2 sets state of EPF2 in next step
new_state_list(:,2) = state_list(:,1);
% Second rule - state of MYC2 sets state of EPFL9 in next step (inhibition)
new_state_list(:,3) = ~state_list(:,1);
% Third rule - only EPF2 and no EPFL9 activates ERECTA
new_state_list( state_list(:,2)&~state_list(:,3),4) = 1;
% Fourth rule - only EPFL9 and no EPF2 inactivates ERECTA
new_state_list(~state_list(:,2)&state_list(:,3),4) = 0;

% Simulate the model
% (Could convert these back to decimal for faster simulation)
cur_state = bool_IC;
bool_X(1,:) = cur_state;
for k = 1:N_gen
    
    % Change to new state based on update rules
    [~,cur_row] = ismember(cur_state, state_list, 'rows');
    cur_state = new_state_list(cur_row,:);
    
    % Store in the trajectory matrix for boolean model
    bool_X(k+1,:) = cur_state;
    
    % Turn off MYC2 when off_gen hit
    if k == off_gen
        cur_state(1) = 0;
    end

end


%%% ODE MODEL

% Define the constants (chosen to match Boolean in this case)
lambda1 = 1; lambda2 = 1; lambda3 = 1; lambda4 = 1;
delta1 = 1; delta2 = 1; delta3 = 1; delta4 = 1;
K1 = 1; K2 = 1; K3 = 1; K4 = 1;

% Define the equations
RHS = @(t,X) [lambda1*(t < off_gen) - delta1*X(1);
             lambda2*X(1)/(K1+X(1)) - delta2*X(2);
             lambda3*K2/(K2+X(1)) - delta3*X(3);
             lambda4*X(2)/(K3+X(2))*K4/(K4+X(3)) - delta4*X(4) ];

% Simulate the equations
ode_options = odeset('AbsTol',1e-6,'RelTol',1e-6);
[T,ode_X] = ode15s( RHS, 0:0.01:N_gen, ode_IC', ode_options);

% Scale all the quantiative model's values to be [0,1]
ode_X = ode_X ./ max(ode_X);

%%% PLOTTING

% Read out useful properties about the problem or the results
N_spec = length(bool_IC);
maxval = max([bool_X(:);ode_X(:)]);

% Plot the trajectories - Boolean in solid lines, ODE in dashed lines
figure('units','Normalized','OuterPosition',[0.2 0.2 0.4 0.6]); hold on;
for k = 1:N_spec
    plot((0:N_gen)-0.05+(k-1)*0.1, bool_X(:,k), 'LineWidth', 4, 'Color', plot_colors(k,:));
    plot(T, ode_X(:,k), '--', 'LineWidth', 5, 'Color', plot_colors(k,:));
end

% Mark the switching off point for MYC2
plot([off_gen off_gen],[0, 1.25*maxval], ':', 'LineWidth', 4, 'Color', [0.6 0.6 0.6]);
annotation('textbox',[off_gen/N_gen-0.2, 0.9*maxval,0.1,0.1], 'String', '(MYC2 On)', 'Color', [0.6 0.6 0.6], 'EdgEColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center');
annotation('textbox',[off_gen/N_gen+0.2, 0.9*maxval,0.1,0.1], 'String', '(MYC2 Off)', 'Color', [0.6 0.6 0.6], 'EdgEColor', 'none', 'FontSize', 22, 'HorizontalAlignment','center');

% Add legend
leg_txt = cell(0);
for k = 1:N_spec
    leg_txt{end+1} = species_names{k};
    leg_txt{end+1} = '';
end
legend(leg_txt, 'FontSize', 24, 'box', 'off', 'Location', 'SouthEast');

% Clean up the plot
xlim([0 N_gen]);
ylim([0 1.05*maxval]);
xlabel('Time');
ylabel('Activity/Concentration');
set(gca,'FontSize',24,'LineWidth',2.5,'Clipping','off');

