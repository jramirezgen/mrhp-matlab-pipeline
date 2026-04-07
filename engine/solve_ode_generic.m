function varargout = solve_ode_generic(cfg, cond_key, sweep_vals)
%SOLVE_ODE_GENERIC Solve one ODE scenario for any organism configuration.
%   result = SOLVE_ODE_GENERIC(cfg, cond_key, sweep_vals)
%   cfg       — full config struct (from organism config file)
%   cond_key  — condition key (e.g. 'MO', 'PHASE1_FUC1', 'CA')
%   sweep_vals — struct with sweep parameter overrides (optional)
%                e.g. struct('substrate_gl',1.0,'target_mM',0.6)

if nargin < 3, sweep_vals = struct(); end

model = cfg.models.(cond_key);
n_sp  = numel(model.species);

% Set initial conditions
x0 = model.base_x0(:);

% Apply sweep overrides OR reference defaults for substrate and target
if isfield(sweep_vals, 'substrate_gl') && isfield(model, 'substrate_idx')
    x0(model.substrate_idx) = sweep_vals.substrate_gl * cfg.solver.substrate_conversion;
elseif isfield(cfg, 'ref_substrate') && isfield(model, 'substrate_idx') && ...
        isfield(cfg.solver, 'substrate_conversion')
    x0(model.substrate_idx) = cfg.ref_substrate * cfg.solver.substrate_conversion;
end
if isfield(sweep_vals, 'target_mM') && isfield(model, 'target_idx')
    x0(model.target_idx) = sweep_vals.target_mM;
elseif isfield(cfg, 'ref_target') && isfield(model, 'target_idx')
    x0(model.target_idx) = cfg.ref_target;
end

% Apply any additional x0 overrides from sweep
if isfield(sweep_vals, 'x0_overrides')
    fns = fieldnames(sweep_vals.x0_overrides);
    for fi = 1:numel(fns)
        idx = find(strcmp(model.species, fns{fi}));
        if ~isempty(idx)
            x0(idx) = sweep_vals.x0_overrides.(fns{fi});
        end
    end
end

% ODE system: dx/dt = S * v(x)
odefun = @(t, x) ode_rhs(t, x, model);

opts = odeset('RelTol', cfg.solver.rtol, 'AbsTol', cfg.solver.atol, ...
              'MaxStep', cfg.solver.max_step, 'NonNegative', 1:n_sp);

try
    [t_sol, y_sol] = ode15s(odefun, cfg.solver.t_eval, x0, opts);
    success = true;
catch ME
    fprintf('  WARN: ODE failed for %s: %s\n', cond_key, ME.message);
    t_sol = cfg.solver.t_eval;
    y_sol = zeros(numel(cfg.solver.t_eval), n_sp);
    success = false;
end

if nargout <= 1
    result.t        = t_sol;
    result.y        = y_sol';        % n_species x n_timepoints
    result.species  = model.species;
    result.condition = cond_key;
    result.success  = success;
    result.x0       = x0;
    % Store sweep values used
    if isfield(sweep_vals, 'substrate_gl')
        result.substrate_gl = sweep_vals.substrate_gl;
    end
    if isfield(sweep_vals, 'target_mM')
        result.target_mM = sweep_vals.target_mM;
    end
    varargout{1} = result;
else
    varargout{1} = t_sol;
    varargout{2} = y_sol;
end
end

function dxdt = ode_rhs(~, x, model)
    x_safe = max(x, 0.0);
    out = model.rate_fn(x_safe, model.params);
    if isstruct(out)
        % Coupled expression: rate_fn returns .v (velocities) + .dxdt_expr
        dxdt_met  = model.S * out.v;
        dxdt = [dxdt_met; out.dxdt_expr];
    else
        % Legacy: rate_fn returns velocity vector directly
        dxdt = model.S * out;
    end
    dxdt(x_safe < 1e-15 & dxdt < 0) = 0;
end
