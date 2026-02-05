% MainModel_batch: runs MainModel_mod.m over a grid of parameter values
%
% Edit struct paramVectors to iterate over; script will run every combination of values and save results to
% 'MainModel_batch_results.mat'.
%

clear;
close all;

% Turn figure visibility off while running the batch (restore after)
oldFigVisible = get(0,'DefaultFigureVisible');
set(0,'DefaultFigureVisible','off');

% -------------------------
% Define parameter vector
% -------------------------
paramVectors = struct();

% EXAMPLE: uncomment / edit to sweep parameters

%%% PARAMETERS WE KEEP FIXED 
paramVectors.r  = 0.05; % discount rate
paramVectors.theta  = [0.971, 0.95];
paramVectors.theta1 = [0.5, 0.6];
paramVectors.yL     = [0.022846, 0.03];
paramVectors.epsiHigh = [0, 0.01];
paramVectors.epsiLow              = 0 ;   % parameter not used for this model

% By default, if no fields set, run a single default case
if isempty(fieldnames(paramVectors))
    disp('No parameter vectors defined in MainModel_batch. Running single default calibration.');
    res = MainModel_main_mod(struct());
    results(1) = res;
else
    fn = fieldnames(paramVectors);
    nFields = numel(fn);

    % Build ndgrid arguments dynamically to create the full grid
    vals = cell(1,nFields);
    for i = 1:nFields
        vals{i} = paramVectors.(fn{i});
    end

    % Use ndgrid to generate full grid
    [gridCells{1:nFields}] = ndgrid(vals{:});
    N = numel(gridCells{1});
    results = repmat(struct(),N,1);

    for k = 1:N
        override = struct();
        for i = 1:nFields
            thisGrid = gridCells{i};
            value = thisGrid(k);
            override.(fn{i}) = value;
        end

        fprintf('Running combination %d/%d ...\n',k,N);
        try
            res = MainModel_main_mod(override);
            results(k) = res;
        catch ME
            % Keep a failure struct for this combination
            results(k).success = false;
            results(k).error = ME;
            results(k).params = override;
            warning('MainModel_batch:runFailed','Run %d failed: %s',k,ME.message);
        end

        % Close any figures generated in this run to free resources
        close all;
    end
end

% Save results to MAT file
save('MainModel_batch_results.mat','results');

% Restore figure visibility
set(0,'DefaultFigureVisible',oldFigVisible);

disp('Batch run finished. Results saved to MainModel_batch_results.mat');