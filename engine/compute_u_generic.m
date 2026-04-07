function u = compute_u_generic(cfg, result, cond_key)
%COMPUTE_U_GENERIC Compute regulatory signals for all genes from config.
%   u = COMPUTE_U_GENERIC(cfg, result, cond_key)
%   cfg      — full config struct
%   result   — ODE result struct from solve_ode_generic
%   cond_key — condition key
%
%   Reads regulatory rules from cfg.regulatory (table of rules per gene).
%   Each rule specifies: form, metabolite_indices, Hill params, basal.

y = result.y;  % n_species x n_time
n_time = size(y, 2);

genes = cfg.genes;
reg   = cfg.regulatory;

% Precompute any derived traces specified in config
derived = struct();
if isfield(cfg, 'derived_signals') && isfield(cfg.derived_signals, cond_key)
    dsigs = cfg.derived_signals.(cond_key);
    fns = fieldnames(dsigs);
    for fi = 1:numel(fns)
        sig_def = dsigs.(fns{fi});
        derived.(fns{fi}) = compute_derived_signal(y, sig_def);
    end
end

% Hill functions
hill_a = @(x, K, n) x.^n ./ (K.^n + x.^n);
hill_r = @(x, K, n) K.^n ./ (K.^n + x.^n);

for gi = 1:numel(genes)
    gene = genes{gi};
    if ~isfield(reg, gene)
        u.(gene) = zeros(1, n_time);
        continue;
    end

    rule = reg.(gene);

    % Get condition-specific rule if available, else use default
    if isfield(rule, cond_key)
        r = rule.(cond_key);
    elseif isfield(rule, 'default')
        r = rule.default;
    else
        r = rule;
    end

    switch r.form
        case 'act'
            % Activation: product of Hill activations over listed signals
            val = ones(1, n_time);
            for si = 1:numel(r.signals)
                sig = get_signal(y, derived, r.signals{si});
                val = val .* hill_a(sig, r.K(si), r.n(si));
            end
            u.(gene) = max(r.basal, val);

        case 'rep'
            % Repression: product of Hill repressions
            val = ones(1, n_time);
            for si = 1:numel(r.signals)
                sig = get_signal(y, derived, r.signals{si});
                val = val .* hill_r(sig, r.K(si), r.n(si));
            end
            u.(gene) = max(r.basal, val);

        case 'act_const'
            % Activation + constitutive: 0.5*Hill + 0.5
            val = zeros(1, n_time);
            for si = 1:numel(r.signals)
                sig = get_signal(y, derived, r.signals{si});
                val = val + hill_a(sig, r.K(si), r.n(si));
            end
            val = val / max(numel(r.signals), 1);
            u.(gene) = max(r.basal, 0.5 * val + 0.5);

        case 'custom'
            % Custom function handle in r.fn
            u.(gene) = r.fn(y, derived, r);

        otherwise
            u.(gene) = r.basal * ones(1, n_time);
    end
end

u.t = result.t(:)';
end

function sig = get_signal(y, derived, signal_spec)
% signal_spec can be:
%   numeric index      → direct species trace: y(idx, :)
%   string 'derived_X' → from derived struct
%   cell {indices}     → sum of species at those indices
    if isnumeric(signal_spec)
        if numel(signal_spec) == 1
            sig = y(signal_spec, :);
        else
            sig = sum(y(signal_spec, :), 1);
        end
    elseif ischar(signal_spec) || isstring(signal_spec)
        sig = derived.(signal_spec);
    elseif iscell(signal_spec)
        sig = sum(y([signal_spec{:}], :), 1);
    else
        sig = zeros(1, size(y, 2));
    end
end

function val = compute_derived_signal(y, sig_def)
% sig_def.type: 'ratio', 'sum', 'diff', 'custom'
% sig_def.indices: species indices involved
    switch sig_def.type
        case 'ratio'
            num = y(sig_def.numerator, :);
            den = sum(y(sig_def.denominator, :), 1);
            den(den < 1e-12) = 1e-12;
            val = num ./ den;
        case 'sum'
            val = sum(y(sig_def.indices, :), 1);
        otherwise
            val = zeros(1, size(y, 2));
    end
end