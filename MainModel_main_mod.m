%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code solves for the stationary equilibrium of the main model
% The parameters are shown in Table 2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = MainModel_main_mod(paramOverrides)

% MainModel_main_mod: Run one instance of the scripts with optional overrides.
%
% paramOverrides - optional struct variable with parameter names to override for this run.
% results - struct with selected outputs and the parameter values used.
%
% Example:
%   res = MainModel_main_mod(struct('theta',0.95,'theta1',0.6));

if nargin < 1
    paramOverrides = struct();
end

% Do not 'clear' here so function arguments remain. Close figures from previous runs.
close all;

% --------------------------
% PARAMETERS AND SETUP 
% --------------------------
LowX                 = 1e-10; %lower bound for integration
TolX                 = 1e-4 ;%tolerance for integration (reduce tolerance to speed up calculations)
gridsize             = 50;  %grid size for the interpolation of value (reduce gridsize to speed up calculations)
xMat                 = linspace(0,1,gridsize) ; %grid of x 
Phi0Mat              = zeros(gridsize,1);
Phi1Mat              = zeros(gridsize,1);
K                    = 7 ; %maximum chain length in markup calculations

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE FUNCTIONS THAT WE NEED IN CALCULATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MainModel_functions ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXOGENOUS PARAMETERS OF THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%% LOAD THE DATA MOMENTS WE KEEP FIXED FOR THE ESTIMATION
MainModel_datamoments;
%%% LOAD A SET OF PARAMETERS
MainModel_paramset ;

% Apply overrides from paramOverrides into the local workspace
overrideFields = fieldnames(paramOverrides);
for i = 1:numel(overrideFields)
    fname = overrideFields{i};
    val = paramOverrides.(fname);
    % assign into local workspace
    eval(sprintf('%s = paramOverrides.%s;', fname, fname));
end

% Start calculations
tic
try
    MainModel_backsolvedemographics ; %identifies demographics parameters
    MainModel_distributions; %solves for distributions
    MainModel_values ; %solve for reservation values
    MainModel_markup ; %solve for markup

    % Capture state after markup but before display (so display can produce figures if desired)
    % display & welfare (these may produce figures & text)
    MainModel_display ; %display the result
    MainModel_price_and_welfare ;

    elapsed = toc;

    % Collect outputs into results struct
    results = struct();

    % Parameters used
    results.params = struct();
    
    % Important params to save 
    Assign = {'theta','theta1','yH','yL','deltaLow','deltaHigh','beta_delta_a','beta_delta_b',...
                   'epsiHigh','beta_epsi_a','beta_epsi_b','M0_over_M','theta0','r','rho','lambda','m','s',...
                   'piH','piL','gamma','m1Eqm','m0Eqm','mul1Eqm','muh0Eqm'};
    for i = 1:numel(Assign)
        name = Assign{i};
        if exist(name,'var')
            results.params.(name) = eval(name);
        end
    end

    % Key equilibrium statistics
    if exist('avg_markup','var'), results.avg_markup = avg_markup; end
    if exist('avg_chain_length','var'), results.avg_chain_length = avg_chain_length; end
    if exist('SizeOfMarkup','var'), results.SizeOfMarkup = SizeOfMarkup; end
    if exist('DistributionOfMarkupMat','var'), results.DistributionOfMarkupMat = DistributionOfMarkupMat; end
    if exist('FrequencyOfChain','var'), results.FrequencyOfChain = FrequencyOfChain; end

    % Welfare & prices if computed by MainModel_price_and_welfare
    if exist('welfare_otc','var'), results.welfare_otc = welfare_otc; end
    if exist('welfare_otc_cust','var'), results.welfare_otc_cust = welfare_otc_cust; end
    if exist('welfare_otc_dealer','var'), results.welfare_otc_dealer = welfare_otc_dealer; end
    if exist('expected_ask','var'), results.expected_ask = expected_ask; end
    if exist('expected_bid','var'), results.expected_bid = expected_bid; end
    if exist('expected_inter','var'), results.expected_inter = expected_inter; end
    if exist('avg_markup_LS','var'), results.avg_markup_LS = avg_markup_LS; end

    results.elapsed_seconds = elapsed;
    results.success = true;

catch ME
    % If an error occurs inside the model calculations, return failure and the error
    elapsed = toc;
    results = struct();
    results.success = false;
    results.error = ME;
    results.elapsed_seconds = elapsed;
    warning('MainModel_main:runError','An error occurred while running the model: %s', ME.message);
end

end