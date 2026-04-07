function phi = extract_phi_generic(cfg, result, cond_key)
%EXTRACT_PHI_GENERIC Extract phi signal from ODE result using config.
%   phi = EXTRACT_PHI_GENERIC(cfg, result, cond_key)
%   Uses cfg.phi_definition to determine how to compute phi.

model = cfg.models.(cond_key);

if isfield(cfg, 'phi_definition') && isfield(cfg.phi_definition, cond_key)
    pdef = cfg.phi_definition.(cond_key);
elseif isfield(cfg, 'phi_definition') && isfield(cfg.phi_definition, 'default')
    pdef = cfg.phi_definition.default;
else
    % Default: 1 - target(t)/target(0)
    pdef = struct('type', 'depletion', 'target_idx', model.target_idx);
end

switch pdef.type
    case 'depletion'
        % phi = 1 - S(t)/S(0)  (like dye decolorization)
        trace = result.y(pdef.target_idx, :);
        s0 = trace(1);
        if s0 <= 0, s0 = result.x0(pdef.target_idx); end
        if s0 <= 0
            phi = zeros(1, size(result.y, 2));
        else
            phi = min(max(1.0 - trace / s0, 0), 1);
        end

    case 'accumulation'
        % phi = P(t)/P_max (like product accumulation)
        trace = result.y(pdef.target_idx, :);
        p_max = pdef.p_max;
        if p_max <= 0, p_max = max(trace); end
        if p_max <= 0
            phi = zeros(1, size(result.y, 2));
        else
            phi = min(max(trace / p_max, 0), 1);
        end

    case 'custom'
        % Custom function handle
        phi = pdef.fn(result, cfg);

    otherwise
        phi = zeros(1, size(result.y, 2));
end
end
