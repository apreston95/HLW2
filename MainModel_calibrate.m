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

    z0 = bounded_to_unbounded(x0, lb, ub);

    obj = @(z) objective_from_unbounded(z, paramNames, lb, ub, targets, calibSpec.display);

    options = optimset('Display', ternary(calibSpec.display,'iter','off'), ...
                       'MaxFunEvals', calibSpec.maxFunEvals, ...
                       'MaxIter', calibSpec.maxIters, ...
                       'TolX', 1e-4, ...
                       'TolFun', 1e-4);

    [zStar, fval, exitflag, output] = fminsearch(obj, z0, options);

    xStar = unbounded_to_bounded(zStar, lb, ub);
    override = build_override(paramNames, xStar);
    bestRun = MainModel_main_mod(override);

    calib = struct();
    calib.paramNames = paramNames;
    calib.bestValues = xStar;
    calib.override = override;
    calib.objective = fval;
    calib.exitflag = exitflag;
    calib.output = output;
    calib.bestRun = bestRun;
    calib.targets = targets;

    if calibSpec.display
        fprintf('\nCalibration complete.\n');
        for i = 1:n
            fprintf('  %s = %.8g\n', paramNames{i}, xStar(i));
        end
        fprintf('  objective = %.8g\n', fval);
    end
end

function loss = objective_from_unbounded(z, paramNames, lb, ub, targets, doDisplay)
    x = unbounded_to_bounded(z, lb, ub);
    override = build_override(paramNames, x);

    run = MainModel_main_mod(override);
    if ~isfield(run,'success') || ~run.success
        loss = 1e9;
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

        scale = max(abs(targetVal), 1e-6);
        err = (modelVal - targetVal) / scale;
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
