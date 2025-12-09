function generateContinuousOutputTables()
% This function prepares the LaTeX tables displayed in the Supplement,
% showing the continuous (non-binned) outputs from the models, relative to
% the wild type value

% Specify how many rows to limit each separate table to
max_rows = 30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load in the particles
load('final_particles.mat','particles');
% Find the best particle
part_Ds = getProperty(particles,'D');
[~,minloc] = min(part_Ds);
best_part = particles{minloc};

% Load in the observed data
[experiments, obs_data] = DunExperiments();
% Prepare the function used to generate model output
f_model = @(theta) runAllExperiments([10.^theta(1:6),theta(7:8)], experiments);
% Define which each dataset corresponds to
data_names = {'cytokinin', 'shoot RMS1', 'root RMS1', 'branching'};
data_filenames = {'CK','RMS1_shoot','RMS1_root','I'};

% Define the base parameter values (values of zero for rates because these
% are log10 of actual value)
base_params = [0, 0, 0, 0, 0, 0, 0.5, 0.5];

% Define the names of the levels
three_level_names = {'Less than WT', 'Similar to WT', 'More than WT'};
five_level_names = {'Much less than WT', 'Less than WT', 'Similar to WT', 'More than WT', 'Much more than WT'};
three_level_short_names = {'Less', 'Similar', 'More'};
five_level_short_names = {'Much less', 'Less', 'Similar', 'More', 'Much more'};

% Run model for the best particle, and also the base (parameter-free) model
best_output = f_model(best_part.theta);
base_output = f_model(base_params);

% Write each table to a text file
for d = 1:length(obs_data)
   
    % Load in experiment names for this dataset (a sheet of the .xlsx file)
    tab_data = readtable('DunData.xlsx','Sheet',data_filenames{d});
    exp_names = tab_data.Description;  
    % Grab out data for this specific set of experiments
    exp_data = obs_data{d};
    best_data = best_output{d};
    base_data = base_output{d};
    % Check using wild type value if there are three or five levels
    if exp_data(1) == 2
        three_levels = true;
        num_levels = 'three';
    elseif exp_data(1) == 3
        three_levels = false;
        num_levels = 'five';
    else
        error('Check wild type value, does not appear to correspond to three or five values');
    end
        
    % Prepare text file
    file_obj = fopen(['continuous_output_',data_filenames{d},'.txt'], 'wt' );
    
    % Each table is split up into multiple pieces if too big, so prepare to
    % scan through all groups of max_rows rows sequentially
    loc = 0;
    table_printed = false;
    looping = true;
    while looping
    
        % Check if this is the final data to print to file
        if size(exp_data,1) - loc <= max_rows
            final_loop = true;
        else
            final_loop = false;
        end
            
        % If this is the final loop, print it as a table, otherwise just a
        % tabular separate from the table environment
        if final_loop
            fprintf(file_obj, '\\begin{table}[!htp]\n');
        else
            fprintf(file_obj, '{\n');    
        end
        
        % Write basic table LaTeX
        fprintf(file_obj, '\\centering\n');
        fprintf(file_obj, '\\begin{tabular}{|c|c|c|c|}\n');
        fprintf(file_obj, '\\hline\n');
        
        % First row is the table 'title' (force capital on first character)
        table_name = [upper(data_names{d}(1)),data_names{d}(2:end)];
        % Append '(Continued)' if we have already printed a table
        if table_printed
            table_name = [table_name, ' (Continued)'];
        else
            table_printed = true;
        end           
        % Print the table name across all rows
        fprintf(file_obj, ['\\multicolumn{4}{|c|}{\\bf ',table_name,'} \\\\ \\hline \n']);
           
        % Next row is the column headers
        fprintf(file_obj, '{\\bf Genotype} & {\\bf Experimental} & {\\bf Param.-free} & {\\bf Best Param.} \\\\ \\hline\n');
    
        % Now loop over current set of rows, successively adding each row
        for r = loc + (1:min([size(exp_data,1) - loc, max_rows]))
       
            %%% Processing of experimental name into LaTeX
            % Read out this experiment's name
            this_exp_name = exp_names{r};
            % Separate out the name into the separators, and texts between
            part_txts = regexp(this_exp_name, '[^/.]+', 'match');
            sep_chars  = regexp(this_exp_name, '[/.]', 'match');
            % Construct the full LaTeX text for experiment name
            exp_txt = '';
            for t = 1:length(part_txts)
                % Add this piece of text, italicising if it is not 'WT'            
                if strcmp(part_txts{t},'WT')
                    exp_txt = [exp_txt, part_txts{t}]; 
                else
                    exp_txt = [exp_txt, '{\\it ', part_txts{t}, '}'];
                end
                % Add the separator character if it exists
                if t <= length(sep_chars)
                    exp_txt = [exp_txt, sep_chars{t}];
                end
            end
        
            %%% Experimental result text
            % If two results are being reported, use short names and separator
            if ~isnan( exp_data(r,2) )
                if three_levels
                    exp_result_txt = [ three_level_short_names{ exp_data(r,1) }, ' | ', three_level_short_names{ exp_data(r,2) } ];
                else
                    exp_result_txt = [ five_level_short_names{ exp_data(r,1) }, ' | ', five_level_short_names{ exp_data(r,2) } ];
                end
            else
                if three_levels
                    exp_result_txt = three_level_names{ exp_data(r,1) };
                else
                    exp_result_txt = five_level_names{ exp_data(r,1) };
                end
            end
        
            %%% Parameter free text
            % Scale model output by wild type value (assumed to come first)
            base_result = base_data(r,:) / base_data(1,1);
            % Convert into text - depends on if two values presented or not
            if ~isnan( exp_data(r,2) )
                base_txt = [num2str(base_result(1), '%.3f'), ' | ', num2str(base_result(2), '%.3f')];
            else
                base_txt = num2str(base_result(1), '%.3f');
            end
        
            %%% Best parameterised model text
            % Scale model output by wild type value (assumed to come first)
            best_result = best_data(r,:) / best_data(1,1);
            % Convert into text - depends on if two values presented or not
            if ~isnan( exp_data(r,2) )
                best_txt = [num2str(best_result(1), '%.3f'), ' | ', num2str(best_result(2), '%.3f')];
            else
                best_txt = num2str(best_result(1), '%.3f');
            end
        
            %%% Adding of row to table
            fprintf(file_obj, [exp_txt, ' & ', exp_result_txt, ' & ', base_txt, ' & ', best_txt, '\\\\\n']);
        
        end
    
        % Finish the LaTeX code for the table itself
        fprintf(file_obj, '\\hline\n');
        fprintf(file_obj, '\\end{tabular}\n');
    
        % If this is the final set of data, add caption and \end{table}
        if final_loop
        
            % Prepare caption text
            caption_txt = ['\\caption{Comparison between experimentally-observed ',data_names{d},' (classified into ',num_levels,' levels)', ...
                           ', and the numerical outputs of the parameter-free and calibrated versions of the Dun {\\it et al.} model relative', ...
                           ' to the wild type (WT) baseline.'];
            % Add additional information to the caption if two measurements shown 
            if any(~isnan(exp_data(:,2)))
                caption_txt = [caption_txt, ' Genotypes that differ between shoot and root are listed shoot/root, and shoot or root grafts', ...
                               ' are indicated by a period. Values separated by | indicate experimental observations (or model predictions)', ...
                               ' for the two grafted scion components.'];
            end
            % Add closing bracket and new line marker
            caption_txt = [caption_txt, '}\n'];
            % Add caption text to the file
            fprintf(file_obj, caption_txt);
    
            % End the table
            fprintf(file_obj, '\\end{table}\n');
           
        % Otherwise, just append a closing bracket and some space
        else
            
            fprintf(file_obj, '}\n\n');
            
        end
    
        % Update location, and terminate loop if this was the final loop
        loc = loc + max_rows;
        if final_loop
            looping = false;
        end
        
    end
    
end
