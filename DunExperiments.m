function [experiments, observations] = DunExperiments()
% This function outputs the experiments, and their observations, contained
% in the Dun et al. data Excel file. Experiments are returned as an array
% (one element per experiment) of structs, where each struct specifies the
% genotypic conditions and the way to convert model output into the
% experimental observation

% Specify effectiveness of genetic knockouts (0 - complete knockout)
RMS1_LEAK = 0;
RMS2_LEAK = 0.2;
RMS3_LEAK = 0;
RMS4_LEAK = 0;
RMS5_LEAK = 0;

% Specify the filename of the datafile
filename = 'DunData.xlsx';

% Get the list of sheet names and the number of experiments
sheets = sheetnames(filename);
N_exp = length(sheets);

% Loop over each sheet in the data
experiments = cell(1,N_exp);
observations = cell(1,N_exp);
for k = 1:N_exp
    
    % Read in the raw data
    data = readtable(filename,'Sheet',sheets(k));
    data_numeric = ensureNumeric(data);
    
    % Trim any completely null rows (all NaN's that come from text rows)
    data_numeric( all( isnan(data_numeric{:,:}), 2),: ) = [];
    
    % Now check data size
    Nrows = size(data_numeric,1);
    
    % Initialise each data
    phi_s = cell(1,Nrows);
    phi_r = cell(1,Nrows);
    results = cell(1,Nrows);
    
    % Read out the experimental configurations and results
    RMS1s1 = data_numeric.RMS1s;
    RMS1s2 = data_numeric.RMS1s2;
    RMS2s1 = data_numeric.RMS2s;
    RMS2s2 = data_numeric.RMS2s2;
    RMS3s1 = data_numeric.RMS3s;
    RMS3s2 = data_numeric.RMS3s2;
    RMS4s1 = data_numeric.RMS4s;
    RMS4s2 = data_numeric.RMS4s2;
    RMS5s1 = data_numeric.RMS5s;
    RMS5s2 = data_numeric.RMS5s2;
    RMS1r1 = data_numeric.RMS1r;
    RMS1r2 = data_numeric.RMS1r2;
    RMS2r1 = data_numeric.RMS2r;
    RMS2r2 = data_numeric.RMS2r2;
    RMS3r1 = data_numeric.RMS3r;
    RMS3r2 = data_numeric.RMS3r2;
    RMS4r1 = data_numeric.RMS4r;
    RMS4r2 = data_numeric.RMS4r2;
    RMS5r1 = data_numeric.RMS5r;
    RMS5r2 = data_numeric.RMS5r2;

    % Replace all "knockout" values in data with their leaky value
    RMS1s1(RMS1s1 == 0) = RMS1_LEAK;
    RMS1s2(RMS1s2 == 0) = RMS1_LEAK;
    RMS1r1(RMS1r1 == 0) = RMS1_LEAK;
    RMS1r2(RMS1r2 == 0) = RMS1_LEAK;
    RMS2s1(RMS2s1 == 0) = RMS2_LEAK;
    RMS2s2(RMS2s2 == 0) = RMS2_LEAK;
    RMS2r1(RMS2r1 == 0) = RMS2_LEAK;
    RMS2r2(RMS2r2 == 0) = RMS2_LEAK;
    RMS3s1(RMS3s1 == 0) = RMS3_LEAK;
    RMS3s2(RMS3s2 == 0) = RMS3_LEAK;
    RMS3r1(RMS3r1 == 0) = RMS3_LEAK;
    RMS3r2(RMS3r2 == 0) = RMS3_LEAK;
    RMS4s1(RMS4s1 == 0) = RMS4_LEAK;
    RMS4s2(RMS4s2 == 0) = RMS4_LEAK;
    RMS4r1(RMS4r1 == 0) = RMS4_LEAK;
    RMS4r2(RMS4r2 == 0) = RMS4_LEAK;
    RMS5s1(RMS5s1 == 0) = RMS5_LEAK;
    RMS5s2(RMS5s2 == 0) = RMS5_LEAK;
    RMS5r1(RMS5r1 == 0) = RMS5_LEAK;
    RMS5r2(RMS5r2 == 0) = RMS5_LEAK;


    %%% SHOOT CONDITIONS AND GENOTYPE
    phi_FS = [RMS2s1, RMS2s2];
    phi_SL = [RMS1s1 .* RMS5s1, RMS1s2 .* RMS5s2];
    phi_I = [RMS3s1 .* RMS4s1, RMS3s2 .* RMS4s2];
    for r = 1:Nrows
        if ~isnan(RMS1s2(r))
            phi_s{r} = [phi_FS(r,:); phi_SL(r,:); phi_I(r,:)];
        else
            phi_s{r} = [phi_FS(r,1); phi_SL(r,1); phi_I(r,1)];
        end
    end
            

    %%% ROOT CONDITIONS AND GENOTYPE
    phi_fs = [RMS2r1, RMS2r2];
    phi_sl = [RMS1r1 .* RMS5r1, RMS1r2 .* RMS5r2];    
    for r = 1:Nrows
        if ~isnan(RMS1r2(r))
            phi_r{r} = [phi_fs(r,:); phi_sl(r,:)];
        else
            phi_r{r} = [phi_fs(r,1); phi_sl(r,1)];
        end
    end   
    
    %%% EXPERIMENTAL RESULT AND SIMULATOR
    results = [data_numeric.result, data_numeric.result2];
    
    % Adjust the [0, 0.5, 1] data to be [1, 2, 3]
    if min(results(:)) == 0
        results = 1 + 2*results;
    end
    
    % Simulator function 
    switch sheets(k)
        
        % Inhibition of branching is fifth row from data
        case 'I'
            for r = 1:Nrows
                f_measure{r} = @(X) X(5,:);
            end
            
        % Cytokinin is sixth row from data
        case 'CK'
            for r = 1:Nrows
                f_measure{r} = @(X) X(6,:);
            end
            
        % RMS1 in shoot is given by FS in shoot and ability to produce RMS1  
        case 'RMS1_shoot'
            for r = 1:Nrows
                if ~isnan(RMS1s2(r))
                    f_measure{r} = @(X) X(1,:) .* [RMS1s1(r), RMS1s2(r)];
                else
                    f_measure{r} = @(X) X(1,:) * RMS1s1(r);
                end
            end
    
        % RMS1 in root is given by FS in root and ability to produce RMS1  
        case 'RMS1_root'
            for r = 1:Nrows
                if ~isnan(RMS1r2(r))
                    f_measure{r} = @(X) X(2,:) .* [RMS1r1(r), RMS1r2(r)];
                else
                    f_measure{r} = @(X) X(2,:) * RMS1r1(r);
                end
            end
        
    end
    
    % Store this experiment in a struct
    experiments{k}.phi_s = phi_s;
    experiments{k}.phi_r = phi_r;
    experiments{k}.results = results;
    experiments{k}.f_measure = f_measure;
    experiments{k}.mutants = [data.Description,data.Description];
    
    % Also store the observations
    observations{k} = results;
    
end


end

function tabdata = ensureNumeric(tabdata)
% This subfunction ensures that all the data inside the provided table is
% numeric

% Read out the size and variable names
[N_rows, N_cols] = size(tabdata);
varnames = tabdata.Properties.VariableNames;

% Loop over each column that does not already contain numeric data
for i = 1:N_cols
    if ~isa(tabdata.(varnames{i}),'double')

        % Create a new numeric vector with empty values to replace with
        vecdata = NaN(N_rows,1);

        % Loop along all column elements and assign them to the vector
        for k = 1:N_rows
            eledata = str2num(tabdata.(varnames{i}){k});
            if ~isempty(eledata)
                vecdata(k) = eledata;
            end
        end

        % Replace table column
        tabdata.(varnames{i}) = vecdata;

    end
end

end
