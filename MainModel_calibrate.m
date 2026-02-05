function calib = MainModel_calibrate(calibSpec)
% MainModel_calibrate: calibrate model parameters to target moments.
%
% Usage:
%   spec = struct();
%   spec.paramNames = {'theta','theta1','yL'};
%   spec.initialValues = [0.97,0.50,0.022846];
%   spec.lowerBounds = [0.85,0.10,0.010];
%   spec.upperBounds = [0.999,0.90,0.060];
%   spec.targets = struct( ...
%       'name',   {'avg_markup','avg_chain_length','expected_inter'}, ...
%       'value',  {2.5,1.35,0.95}, ...
%       'weight', {1.0,1.0,1.0});
%   spec.display = true;
%   calib = MainModel_calibrate(spec);
%
% Required calibSpec fields:
%   - paramNames (cellstr)
%   - initialValues (numeric vector)
%   - lowerBounds (numeric vector)
%   - upperBounds (numeric vector)
%   - targets (struct array with fields: name, value, weight)
%
% Optional fields:
%   - display (logical, default=true)
%   - maxFunEvals (default=300)
%   - maxIters (default=200)
%   - nStarts (default=3): number of optimization restarts
%   - feasibleSearchDraws (default=40): random points used if x0 fails
%
% Returns struct calib with best params, objective value, optimizer output,
% and model results at optimum.

    validateattributes(calibSpec, {'struct'}, {'nonempty'});

    required = {'paramNames','initialValues','lowerBounds','upperBounds','targets'};
    for i = 1:numel(required)
        if ~isfield(calibSpec, required{i})
            error('MainModel_calibrate:missingField', 'Missing calibSpec.%s', required{i});
        end
    end

    paramNames = calibSpec.paramNames;
    x0 = calibSpec.initialValues(:)';
    lb = calibSpec.lowerBounds(:)';
    ub = calibSpec.upperBounds(:)';
    targets = calibSpec.targets;

    n = numel(paramNames);
    if numel(x0) ~= n || numel(lb) ~= n || numel(ub) ~= n
        error('MainModel_calibrate:dimensionMismatch', ...
            'paramNames, initialValues, lowerBounds and upperBounds must have same length.');
    end

    if any(lb >= ub)
        error('MainModel_calibrate:invalidBounds','Each lower bound must be < upper bound.');
    end

    if ~isstruct(targets) || isempty(targets)
        error('MainModel_calibrate:invalidTargets','targets must be a non-empty struct array.');
    end

    if ~isfield(calibSpec,'display'), calibSpec.display = true; end
    if ~isfield(calibSpec,'maxFunEvals'), calibSpec.maxFunEvals = 300; end
    if ~isfield(calibSpec,'maxIters'), calibSpec.maxIters = 200; end
    if ~isfield(calibSpec,'nStarts'), calibSpec.nStarts = 3; end
    if ~isfield(calibSpec,'feasibleSearchDraws'), calibSpec.feasibleSearchDraws = 40; end

    % Ensure we start from a point where the model can actually run.
    xStart = find_feasible_start(x0, lb, ub, paramNames, calibSpec.feasibleSearchDraws, calibSpec.display);

    options = optimset('Display', ternary(calibSpec.display,'iter','off'), ...
                       'MaxFunEvals', calibSpec.maxFunEvals, ...
                       'MaxIter', calibSpec.maxIters, ...
                       'TolX', 1e-4, ...
                       'TolFun', 1e-4);

    best = struct('fval', inf, 'z', [], 'x', [], 'exitflag', [], 'output', []);

    for s = 1:calibSpec.nStarts
        if s == 1
            xSeed = xStart;
        else
            jitter = (ub-lb) .* (0.10 .* (2*rand(1,n)-1));
            xSeed = min(max(xStart + jitter, lb), ub);
        end

        z0 = bounded_to_unbounded(xSeed, lb, ub);
        obj = @(z) objective_from_unbounded(z, paramNames, lb, ub, targets, calibSpec.display);
        [zStar, fval, exitflag, output] = fminsearch(obj, z0, options);

        if fval < best.fval
            best.fval = fval;
            best.z = zStar;
            best.x = unbounded_to_bounded(zStar, lb, ub);
            best.exitflag = exitflag;
            best.output = output;
        end
    end

    xStar = best.x;
    override = build_override(paramNames, xStar);
    override = apply_dependent_overrides(override);
    bestRun = MainModel_main_mod(override);

    calib = struct();
    calib.paramNames = paramNames;
    calib.bestValues = xStar;
    calib.override = override;
    calib.objective = best.fval;
    calib.exitflag = best.exitflag;
    calib.output = best.output;
    calib.bestRun = bestRun;
    calib.targets = targets;

    if calibSpec.display
        fprintf('\nCalibration complete.\n');
        for i = 1:n
            fprintf('  %s = %.8g\n', paramNames{i}, xStar(i));
        end
        fprintf('  objective = %.8g\n', best.fval);
    end
end

function xStart = find_feasible_start(x0, lb, ub, paramNames, nDraws, doDisplay)
    [isFeasible, ~] = check_feasible(x0, paramNames);
    if isFeasible
        xStart = x0;
        return;
    end

    if doDisplay
        fprintf('Initial values are infeasible. Searching for a feasible start...\n');
    end

    xStart = [];
    for k = 1:nDraws
        u = rand(1, numel(x0));
        cand = lb + (ub-lb).*u;
        [ok, ~] = check_feasible(cand, paramNames);
        if ok
            xStart = cand;
            break;
        end
    end

    if isempty(xStart)
        error('MainModel_calibrate:noFeasibleStart', ...
            ['No feasible start found in bounds. ', ...
             'Model runs are failing for tested points; expand bounds or adjust parameters.']);
    end
end

function [ok, run] = check_feasible(x, paramNames)
    override = build_override(paramNames, x);
    override = apply_dependent_overrides(override);
    run = MainModel_main_mod(override);
    ok = isfield(run,'success') && run.success;
end

function loss = objective_from_unbounded(z, paramNames, lb, ub, targets, doDisplay)
    x = unbounded_to_bounded(z, lb, ub);
    override = build_override(paramNames, x);
    override = apply_dependent_overrides(override);

    run = MainModel_main_mod(override);
    if ~isfield(run,'success') || ~run.success
        % Non-constant penalty to help optimizer move off bad regions.
        mid = 0.5*(lb+ub);
        scale = max(ub-lb, 1e-6);
        loss = 1e9 + 1e3*sum(((x-mid)./scale).^2);
        return;
    end

    loss = 0;
    for i = 1:numel(targets)
        name = targets(i).name;
        targetVal = targets(i).value;
        weight = targets(i).weight;

        if ~isfield(run, name)
            loss = loss + 1e6;
            continue;
        end

        modelVal = run.(name);
        if ~isfinite(modelVal)
            loss = loss + 1e6;
            continue;
        end

        s = max(abs(targetVal), 1e-6);
        err = (modelVal - targetVal) / s;
        loss = loss + weight * err^2;
    end

    if doDisplay
        fprintf('loss=%g | ', loss);
        for j = 1:numel(paramNames)
            fprintf('%s=%g ', paramNames{j}, x(j));
        end
        fprintf('\n');
    end
end

function z = bounded_to_unbounded(x, lb, ub)
    xn = (x - lb) ./ (ub - lb);
    xn = min(max(xn, 1e-8), 1-1e-8);
    z = log(xn ./ (1 - xn));
end

function x = unbounded_to_bounded(z, lb, ub)
    xn = 1 ./ (1 + exp(-z));
    x = lb + (ub - lb) .* xn;
end

function override = build_override(paramNames, x)
    override = struct();
    for i = 1:numel(paramNames)
        override.(paramNames{i}) = x(i);
    end
end

function override = apply_dependent_overrides(override)
    % Keep internally linked parameters coherent during optimization.
    if isfield(override,'theta1') && ~isfield(override,'theta0')
        override.theta0 = 1 - override.theta1;
    end
    if isfield(override,'r') && ~isfield(override,'yH')
        override.yH = override.r;
    end
end

function y = ternary(cond, a, b)
    if cond
        y = a;
    else
        y = b;
    end
end
