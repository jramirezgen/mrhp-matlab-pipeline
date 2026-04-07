function D = phenotype_reference_generic(cfg, cond_key, sweep_vals, t)
%PHENOTYPE_REFERENCE_GENERIC Compute phenomenological reference curve.
%   D = PHENOTYPE_REFERENCE_GENERIC(cfg, cond_key, sweep_vals, t)
%   Uses cfg.phenotype_ref to generate the reference curve.

if nargin < 3, sweep_vals = struct(); end
if nargin < 4, t = cfg.solver.t_eval; end

pref = cfg.phenotype_ref.(cond_key);

switch pref.type
    case 'hill_decolorization'
        % V4-style: D = Dmax * (1 - exp(-Vmax * h_SI * t_eff))
        S = pref.default_substrate;
        I = pref.default_target;
        if isfield(sweep_vals, 'substrate_gl'), S = sweep_vals.substrate_gl; end
        if isfield(sweep_vals, 'target_mM'),    I = sweep_vals.target_mM; end
        h_SI = (S^pref.n / (pref.K_yc^pref.n + S^pref.n)) * (pref.K_dic / (pref.K_dic + I));
        t_eff = max(t - pref.lag, 0);
        D = pref.Dmax * (1 - exp(-pref.Vmax * h_SI * t_eff));
        D = min(max(D, 0), 100);

    case 'sigmoidal_accumulation'
        % P(t) = Pmax * (1 - exp(-k * (t - lag)^n))
        t_eff = max(t - pref.lag, 0);
        D = pref.Pmax * (1 - exp(-pref.k * t_eff .^ pref.n));
        D = min(max(D, 0), pref.Pmax);

    case 'exponential_depletion'
        % S(t) = S0 * exp(-mu * t) → D = (1 - S/S0)*100
        D = 100 * (1 - exp(-pref.mu * t));
        D = min(max(D, 0), 100);

    case 'data_interpolation'
        % Interpolate from experimental data
        D = interp1(pref.t_data, pref.y_data, t, 'pchip', pref.y_data(end));

    case 'custom'
        D = pref.fn(t, sweep_vals, cfg);

    otherwise
        D = zeros(size(t));
end
end
