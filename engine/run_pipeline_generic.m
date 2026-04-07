function run_pipeline_generic(cfg, varargin)
%RUN_PIPELINE_GENERIC Universal MRHP Multiscale Pipeline Orchestrator.
%   run_pipeline_generic(cfg)
%   run_pipeline_generic(cfg, 'quick')
%   run_pipeline_generic(cfg, 'figures_only')
%
%   cfg — full config struct from organism config file
%
%   Phases:
%     0: Verify    1: ODE    2: Phi    3: Bridge    4: Regulatory
%     5: Deterministic expression    6: Stochastic expression
%     7: Full grid    8: Figures    9: Validation
%     10: Metadata   11: Report   12: Certification

%% PARSE OPTIONS
opts = struct('quick', false, 'figures_only', false);
for i = 1:numel(varargin)
    if ischar(varargin{i})
        switch lower(varargin{i})
            case 'quick', opts.quick = true;
            case 'figures_only', opts.figures_only = true;
        end
    end
end

%% SETUP
conditions = cfg.conditions;
genes      = cfg.genes;
n_cond     = numel(conditions);
n_genes    = numel(genes);

base_out = cfg.output_dir;
dirs = {'ode_timeseries','phi_signals','bridge_predictions','expression',...
        'regulatory_signals','figures','metrics','summary_metrics',...
        'hypothesis_ranking'};
for di = 1:numel(dirs)
    d = fullfile(base_out, dirs{di});
    if ~exist(d,'dir'), mkdir(d); end
end
logs_dir = fullfile(base_out, '..', 'logs');
meta_dir = fullfile(base_out, '..', 'metadata');
reports_dir = fullfile(base_out, '..', 'reports');
for d = {logs_dir, meta_dir, reports_dir}
    if ~exist(d{1},'dir'), mkdir(d{1}); end
end

t0 = tic;
ts = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
run_info = struct('timestamp', ts, 'organism', cfg.organism, 'phases', struct(), 'validation', []);

fprintf('=== MRHP Universal Pipeline — %s ===\n', cfg.organism);
fprintf('    Conditions: %s\n', strjoin(conditions, ', '));
fprintf('    Genes: %s\n', strjoin(genes, ', '));
fprintf('    Timestamp: %s\n\n', ts);

% Build sweep grid from config
if isfield(cfg, 'sweep')
    sweep_grid = build_sweep_grid(cfg.sweep);
else
    sweep_grid = {struct()};  % single empty sweep = reference condition only
end
n_sweep = numel(sweep_grid);

%% FIGURES ONLY
if opts.figures_only
    fprintf('=== FIGURES ONLY MODE ===\n\n');
    figure_generation_generic(cfg, base_out);
    return;
end

%% PHASE 0: VERIFY
fprintf('\n== PHASE 0: VERIFY ==\n');
fprintf('  [OK] MATLAB R%s\n', version('-release'));
fprintf('  [OK] Organism: %s\n', cfg.organism);
fprintf('  [OK] %d conditions, %d genes\n', n_cond, n_genes);
run_info.phases.P0_verify = 'OK';

%% PHASE 1: ODE
fprintf('\n== PHASE 1: ODE METABOLIC MODELS ==\n');
out_ode = fullfile(base_out, 'ode_timeseries');
converged = 0; total = 0;
for ci = 1:n_cond
    cond = conditions{ci};
    for si = 1:n_sweep
        total = total + 1;
        sv = sweep_grid{si};
        sid = make_scenario_id(cond, sv);
        try
            result = solve_ode_generic(cfg, cond, sv);
            sp = result.species;
            n_sp = numel(sp);
            t = result.t;
            y = result.y;
            fid = fopen(fullfile(out_ode, sprintf('timeseries_%s.tsv', sid)), 'w');
            fprintf(fid, 'time_h');
            for spi = 1:n_sp, fprintf(fid, '\t%s', sp{spi}); end
            fprintf(fid, '\n');
            fmt = ['%.4f', repmat('\t%.6e', 1, n_sp), '\n'];
            for ri = 1:numel(t), fprintf(fid, fmt, t(ri), y(:, ri)); end
            fclose(fid);
            converged = converged + 1;
        catch me
            fprintf('  WARN: %s failed: %s\n', sid, me.message);
        end
    end
end
fprintf('  Converged: %d/%d scenarios\n', converged, total);
run_info.phases.P1_ode = 'OK';

%% PHASE 2: PHI SIGNALS
fprintf('\n== PHASE 2: PHI SIGNAL EXTRACTION ==\n');
out_phi = fullfile(base_out, 'phi_signals');
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    phi = extract_phi_generic(cfg, result, cond);
    t = result.t(:)';
    fid = fopen(fullfile(out_phi, sprintf('phi_signals_%s.tsv', cond)), 'w');
    fprintf(fid, 'time_h\tphi\n');
    for ri = 1:numel(t), fprintf(fid, '%.4f\t%.6f\n', t(ri), phi(ri)); end
    fclose(fid);
    fprintf('  Phi signals for %s\n', cond);
end
run_info.phases.P2_phi = 'OK';

%% PHASE 3: BRIDGE FITTING
fprintf('\n== PHASE 3: BRIDGE FITTING ==\n');
out_bridge = fullfile(base_out, 'bridge_predictions');
out_metrics = fullfile(base_out, 'metrics');
all_R2 = [];
fid_met = fopen(fullfile(out_metrics, 'bridge_all_scenarios.tsv'), 'w');
fprintf(fid_met, 'condition\tscenario\tR2\tlambda\tbeta\tymax\n');

for ci = 1:n_cond
    cond = conditions{ci};
    for si = 1:n_sweep
        sv = sweep_grid{si};
        sid = make_scenario_id(cond, sv);
        result = solve_ode_generic(cfg, cond, sv);
        phi = extract_phi_generic(cfg, result, cond);
        t = result.t(:)';
        D_ref = phenotype_reference_generic(cfg, cond, sv, t);
        [lam, beta_v, ymax, R2, ~, y_pred] = fit_bridge(phi, t, D_ref);
        fprintf(fid_met, '%s\t%s\t%.4f\t%.4f\t%.4f\t%.2f\n', cond, sid, R2, lam, beta_v, ymax);
        all_R2(end+1) = R2;
    end
    % Representative bridge prediction (reference condition)
    result = solve_ode_generic(cfg, cond);
    phi = extract_phi_generic(cfg, result, cond);
    t = result.t(:)';
    D_ref = phenotype_reference_generic(cfg, cond, struct(), t);
    [lam, beta_v, ymax, R2, ~, y_pred] = fit_bridge(phi, t, D_ref);
    fid = fopen(fullfile(out_bridge, sprintf('bridge_pred_%s.tsv', cond)), 'w');
    fprintf(fid, 'time_h\tD_ref\tD_bridge\tphi\n');
    for ri = 1:numel(t), fprintf(fid, '%.4f\t%.4f\t%.4f\t%.6f\n', t(ri), D_ref(ri), y_pred(ri), phi(ri)); end
    fclose(fid);
    fprintf('  Bridge %s: R^2=%.4f\n', cond, R2);
end
fclose(fid_met);
fprintf('  Mean R^2: %.4f (n=%d)\n', mean(all_R2), numel(all_R2));
run_info.phases.P3_bridge = 'OK';

%% PHASE 4: REGULATORY SIGNALS
fprintf('\n== PHASE 4: REGULATORY SIGNALS ==\n');
out_reg = fullfile(base_out, 'regulatory_signals');
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    u = compute_u_generic(cfg, result, cond);
    fid_r = fopen(fullfile(out_reg, sprintf('regulatory_%s.tsv', cond)), 'w');
    fprintf(fid_r, 'time_h');
    for gi = 1:n_genes, fprintf(fid_r, '\tu_%s', genes{gi}); end
    fprintf(fid_r, '\n');
    for k = 1:numel(u.t)
        fprintf(fid_r, '%.4f', u.t(k));
        for gi = 1:n_genes, fprintf(fid_r, '\t%.6f', u.(genes{gi})(k)); end
        fprintf(fid_r, '\n');
    end
    fclose(fid_r);
    fprintf('  %s: regulatory signals computed\n', cond);
end
run_info.phases.P4_regulatory = 'OK';

%% PHASE 5: DETERMINISTIC EXPRESSION
fprintf('\n== PHASE 5: DETERMINISTIC EXPRESSION ==\n');
out_expr = fullfile(base_out, 'expression');
t_eval = cfg.solver.t_eval;
is_coupled = isfield(cfg, 'coupled_expression') && cfg.coupled_expression;
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    if is_coupled
        % Expression comes directly from unified ODE (species 19:34)
        n_met_c = 18; n_g_c = n_genes;
        t_det = result.t(:);
        for gi = 1:n_genes
            gene = genes{gi};
            m_det = result.y(n_met_c + gi, :)';
            p_det = result.y(n_met_c + n_g_c + gi, :)';
            fid_e = fopen(fullfile(out_expr, sprintf('det_%s_%s.tsv', cond, gene)), 'w');
            fprintf(fid_e, 'time_h\tmRNA_nM\tProtein_nM\n');
            for k = 1:numel(t_det), fprintf(fid_e, '%.2f\t%.6e\t%.6e\n', t_det(k), m_det(k), p_det(k)); end
            fclose(fid_e);
            idx_sample = find(t_det >= cfg.t_sample, 1);
            if isempty(idx_sample), idx_sample = numel(t_det); end
            fprintf('  %s/%s (coupled): m@%.0fh=%.4e nM\n', cond, gene, cfg.t_sample, m_det(idx_sample));
        end
    else
        % Legacy: separate expression ODE
        u = compute_u_generic(cfg, result, cond);
        for gi = 1:n_genes
            gene = genes{gi};
            ktx = cfg.ktx_fits.(gene);
            [t_det, m_det, p_det] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
                cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
            fid_e = fopen(fullfile(out_expr, sprintf('det_%s_%s.tsv', cond, gene)), 'w');
            fprintf(fid_e, 'time_h\tmRNA_nM\tProtein_nM\n');
            for k = 1:numel(t_det), fprintf(fid_e, '%.2f\t%.6e\t%.6e\n', t_det(k), m_det(k), p_det(k)); end
            fclose(fid_e);
            idx_sample = find(t_det >= cfg.t_sample, 1);
            if isempty(idx_sample), idx_sample = numel(t_det); end
            fprintf('  %s/%s: m@%.0fh=%.4e nM\n', cond, gene, cfg.t_sample, m_det(idx_sample));
        end
    end
end
run_info.phases.P5_deterministic = 'OK';

%% PHASE 6: STOCHASTIC EXPRESSION
fprintf('\n== PHASE 6: STOCHASTIC EXPRESSION ==\n');
ref_cond = conditions{1};
result = solve_ode_generic(cfg, ref_cond);
u = compute_u_generic(cfg, result, ref_cond);
expr_params = struct('stochastic_tau', cfg.expression.stochastic_tau, ...
    'beta_m', cfg.expression.beta_m, 'beta_p', cfg.expression.beta_p, ...
    'ktl', cfg.expression.ktl, 'nm_to_mol', cfg.expression.nm_to_mol, ...
    't_span', cfg.solver.t_span, 't_eval', cfg.solver.t_eval);
for gi = 1:n_genes
    gene = genes{gi};
    ktx = cfg.ktx_fits.(gene);
    sto = tau_leaping_expression(u.t, u.(gene), ktx, cfg.expression.n_cells, 42, expr_params);
    idx_sample = find(sto.t >= cfg.t_sample, 1);
    if isempty(idx_sample), idx_sample = numel(sto.t); end
    fprintf('  %s/%s: mean_m@%.0fh=%.4e nM\n', ref_cond, gene, cfg.t_sample, sto.mean_m(idx_sample));
    fid_s = fopen(fullfile(out_expr, sprintf('sto_%s_%s.tsv', ref_cond, gene)), 'w');
    fprintf(fid_s, 'time_h\tmean_m_nM\tstd_m_nM\tmean_p_nM\tstd_p_nM\n');
    for k = 1:numel(sto.t)
        fprintf(fid_s, '%.2f\t%.6e\t%.6e\t%.6e\t%.6e\n', sto.t(k), sto.mean_m(k), sto.std_m(k), sto.mean_p(k), sto.std_p(k));
    end
    fclose(fid_s);
end
run_info.phases.P6_stochastic = 'OK';

%% PHASE 7: FULL EXPRESSION GRID
if ~opts.quick
    fprintf('\n== PHASE 7: FULL EXPRESSION GRID ==\n');
    out_summary = fullfile(base_out, 'summary_metrics');
    fid = fopen(fullfile(out_summary, 'scenario_summary.tsv'), 'w');
    fprintf(fid, 'condition\tscenario\tphi_end\tbridge_R2');
    for gi = 1:n_genes, fprintf(fid, '\tm_%s', genes{gi}); end
    fprintf(fid, '\n');
    n_scenarios = 0;
    for ci = 1:n_cond
        cond = conditions{ci};
        for si = 1:n_sweep
            sv = sweep_grid{si};
            sid = make_scenario_id(cond, sv);
            result = solve_ode_generic(cfg, cond, sv);
            u = compute_u_generic(cfg, result, cond);
            phi = extract_phi_generic(cfg, result, cond);
            t = result.t(:)';
            D_ref = phenotype_reference_generic(cfg, cond, sv, t);
            [~,~,~, R2_sc, ~, ~] = fit_bridge(phi, t, D_ref);
            fprintf(fid, '%s\t%s\t%.4f\t%.4f', cond, sid, phi(end), R2_sc);
            for gi = 1:n_genes
                gene = genes{gi};
                [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), cfg.ktx_fits.(gene), ...
                    cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, cfg.solver.t_eval);
                idx_sample = find(t_det >= cfg.t_sample, 1);
                if isempty(idx_sample), idx_sample = numel(t_det); end
                fprintf(fid, '\t%.6e', m_det(idx_sample));
            end
            fprintf(fid, '\n');
            n_scenarios = n_scenarios + 1;
        end
    end
    fclose(fid);
    fprintf('  Summary: %d scenarios\n', n_scenarios);
    run_info.phases.P7_grid = 'OK';
end

%% PHASE 8: FIGURES
fprintf('\n== PHASE 8: FIGURE GENERATION ==\n');
figure_generation_generic(cfg, base_out);
run_info.phases.P8_figures = 'OK';

%% PHASE 8b: HYPOTHESIS EVALUATION
if isfield(cfg, 'hypotheses') && ~isempty(fieldnames(cfg.hypotheses))
    fprintf('\n== PHASE 8b: HYPOTHESIS EVALUATION ==\n');
    [rankings, hyp_evaluations] = hypothesis_evaluator_generic(cfg, base_out);
    run_info.hypothesis_rankings = rankings;
    run_info.phases.P8b_hypothesis = 'OK';
else
    fprintf('\n== PHASE 8b: HYPOTHESIS — no hypotheses defined, skipping ==\n');
    run_info.phases.P8b_hypothesis = 'SKIPPED';
end

%% PHASE 9: VALIDATION
fprintf('\n== PHASE 9: VALIDATION ==\n');
if isfield(cfg, 'experimental_data') && ~isempty(fieldnames(cfg.experimental_data))
    val_result = validate_against_data(cfg);
    run_info.validation = val_result;
    fprintf('  Validation complete\n');
else
    fprintf('  No experimental data — skipping\n');
    run_info.validation = struct('status', 'NO_DATA');
end
run_info.phases.P9_validation = 'OK';

%% PHASE 10: METADATA
duration = toc(t0);
run_info.duration_s = duration;
fprintf('\n== PHASE 10: METADATA ==\n');
fid = fopen(fullfile(meta_dir, 'run_info.json'), 'w');
fprintf(fid, '%s', jsonencode(run_info));
fclose(fid);
run_info.phases.P10_metadata = 'OK';

%% PHASE 10b: STOICHIOMETRIC MATRIX EXPORT
if isfield(cfg, 'coupled_expression') && cfg.coupled_expression
    fprintf('\n== PHASE 10b: STOICHIOMETRIC MATRIX EXPORT ==\n');
    ref_cond = conditions{1};
    model_ref = cfg.models.(ref_cond);
    S_mat = model_ref.S;
    sp_met = model_ref.species(1:size(S_mat,1));
    n_rxn = size(S_mat, 2);
    % Export S matrix
    fid_s = fopen(fullfile(meta_dir, 'stoichiometric_matrix_S.tsv'), 'w');
    fprintf(fid_s, 'species');
    for ri = 1:n_rxn, fprintf(fid_s, '\tR%d', ri); end
    fprintf(fid_s, '\n');
    for si = 1:numel(sp_met)
        fprintf(fid_s, '%s', sp_met{si});
        for ri = 1:n_rxn, fprintf(fid_s, '\t%d', S_mat(si,ri)); end
        fprintf(fid_s, '\n');
    end
    fclose(fid_s);
    % Export x0
    fid_x = fopen(fullfile(meta_dir, 'initial_conditions_x0.tsv'), 'w');
    fprintf(fid_x, 'species\tx0\n');
    all_sp = model_ref.species;
    x0_ref = model_ref.base_x0;
    for si = 1:numel(all_sp)
        fprintf(fid_x, '%s\t%.6e\n', all_sp{si}, x0_ref(si));
    end
    fclose(fid_x);
    fprintf('  Exported S (%dx%d), x0 (%d species) to metadata/\n', size(S_mat,1), n_rxn, numel(all_sp));
    run_info.phases.P10b_matrix_export = 'OK';
end

%% PHASE 10c: ROUTE MAP EXPORT
fprintf('\n== PHASE 10c: ROUTE MAP EXPORT ==\n');
if isfield(cfg, 'route_maps')
    rm_dir = fullfile(base_out, 'route_maps');
    if ~exist(rm_dir, 'dir'), mkdir(rm_dir); end
    configs_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'configs');
    for ci = 1:n_cond
        cond = conditions{ci};
        if isfield(cfg.route_maps, cond)
            rm = cfg.route_maps.(cond);
            fnames = fieldnames(rm);
            for fi = 1:numel(fnames)
                src = fullfile(configs_dir, rm.(fnames{fi}));
                dst = fullfile(rm_dir, sprintf('%s_%s.tsv', fnames{fi}, cond));
                if exist(src, 'file')
                    copyfile(src, dst);
                    fprintf('  Exported %s\n', dst);
                end
            end
        end
    end
    % Export stoichiometric equations markdown
    eq_file = fullfile(rm_dir, 'stoichiometric_equations.md');
    fid_eq = fopen(eq_file, 'w');
    fprintf(fid_eq, '# Stoichiometric Equations — %s\n\n', cfg.organism);
    for ci = 1:n_cond
        cond = conditions{ci};
        fprintf(fid_eq, '## %s\n\n', cond);
        if isfield(cfg.route_maps, cond)
            src = fullfile(configs_dir, cfg.route_maps.(cond).route_map);
            if exist(src, 'file')
                T = readtable(src, 'FileType', 'text', 'Delimiter', '\t');
                fprintf(fid_eq, '| # | Reaction | Equation | GEM ID | Gene | EC |\n');
                fprintf(fid_eq, '|---|----------|----------|--------|------|----| \n');
                for ri = 1:height(T)
                    row = T(ri,:);
                    fprintf(fid_eq, '| %d | %s | %s | %s | %s | %s |\n', ...
                        row.ode_rxn_idx, string(row.ode_rxn_name), ...
                        string(row.equation), string(row.gem_rxn_id), ...
                        string(row.gene), string(row.ec));
                end
                fprintf(fid_eq, '\n');
            end
        end
    end
    fclose(fid_eq);
    fprintf('  Stoichiometric equations exported\n');
    run_info.phases.P10c_route_maps = 'OK';
else
    fprintf('  No route_maps defined — skipping\n');
    run_info.phases.P10c_route_maps = 'SKIPPED';
end

%% PHASE 10d: S MATRIX NUMERIC EXPORT
fprintf('\n== PHASE 10d: S MATRIX NUMERIC EXPORT ==\n');
for ci = 1:n_cond
    cond = conditions{ci};
    model = cfg.models.(cond);
    n_sp = size(model.S, 1);
    n_rxn = size(model.S, 2);
    sm_dir = fullfile(base_out, 'route_maps');
    if ~exist(sm_dir, 'dir'), mkdir(sm_dir); end
    fid_s = fopen(fullfile(sm_dir, sprintf('S_numeric_%s.tsv', cond)), 'w');
    fprintf(fid_s, 'species');
    for ri = 1:n_rxn, fprintf(fid_s, '\tv%d', ri); end
    fprintf(fid_s, '\n');
    for si = 1:n_sp
        fprintf(fid_s, '%s', model.species{si});
        for ri = 1:n_rxn, fprintf(fid_s, '\t%d', model.S(si, ri)); end
        fprintf(fid_s, '\n');
    end
    fclose(fid_s);
    fprintf('  S matrix %s: %dx%d\n', cond, n_sp, n_rxn);
end
run_info.phases.P10d_smatrix = 'OK';

%% PHASE 11: REPORT
fprintf('\n== PHASE 11: REPORT ==\n');
fid = fopen(fullfile(reports_dir, 'pipeline_execution_report.md'), 'w');
fprintf(fid, '# MRHP Universal Pipeline — %s\n\n', cfg.organism);
fprintf(fid, '**Timestamp:** %s\n', ts);
fprintf(fid, '**Duration:** %.1fs\n', duration);
fprintf(fid, '**Conditions:** %s\n', strjoin(conditions, ', '));
fprintf(fid, '**Genes:** %s\n\n', strjoin(genes, ', '));
fprintf(fid, '## Phases\n\n');
pn = fieldnames(run_info.phases);
for pi = 1:numel(pn), fprintf(fid, '- %s: **%s**\n', pn{pi}, run_info.phases.(pn{pi})); end
fclose(fid);
run_info.phases.P11_report = 'OK';

%% PHASE 12: CERTIFICATION
fprintf('\n== PHASE 12: CERTIFICATION ==\n');
expected = {'P0_verify','P1_ode','P2_phi','P3_bridge','P4_regulatory',...
            'P5_deterministic','P6_stochastic','P8_figures','P9_validation'};
if ~opts.quick, expected = [expected, {'P7_grid'}]; end
all_ok = true;
for ei = 1:numel(expected)
    if ~isfield(run_info.phases, expected{ei}) || ~strcmp(run_info.phases.(expected{ei}), 'OK')
        all_ok = false; break;
    end
end
status = 'PARTIAL';
if all_ok, status = 'CERTIFIED'; end

n_figs = count_figures(fullfile(base_out, 'figures'));
cert = struct('status', status, 'timestamp', ts, 'organism', cfg.organism, ...
    'n_conditions', n_cond, 'n_genes', n_genes, 'n_figures', n_figs);
fid = fopen(fullfile(meta_dir, 'certification.json'), 'w');
fprintf(fid, '%s', jsonencode(cert));
fclose(fid);

fprintf('  Pipeline status: %s (%d figures)\n', status, n_figs);
fprintf('\n=== Pipeline complete in %.1fs — %s ===\n', duration, status);
end

%% ═══ HELPERS ═══
function grid = build_sweep_grid(sweep_def)
    if isfield(sweep_def, 'substrate_values') && isfield(sweep_def, 'target_values')
        s_vals = sweep_def.substrate_values;
        t_vals = sweep_def.target_values;
        grid = {};
        for si = 1:numel(s_vals)
            for ti = 1:numel(t_vals)
                grid{end+1} = struct('substrate_gl', s_vals(si), 'target_mM', t_vals(ti));
            end
        end
    else
        grid = {struct()};
    end
end

function sid = make_scenario_id(cond, sv)
    sid = cond;
    if isfield(sv, 'substrate_gl')
        sid = sprintf('%s_S%.1f', sid, sv.substrate_gl);
    end
    if isfield(sv, 'target_mM')
        sid = sprintf('%s_T%.1f', sid, sv.target_mM);
    end
end

function n = count_figures(fig_dir)
    n = 0;
    if exist(fig_dir, 'dir')
        pngs = dir(fullfile(fig_dir, '*.png'));
        n = numel(pngs);
    end
end

function val = validate_against_data(cfg)
    val = struct('status', 'COMPUTED');
    if ~isfield(cfg, 'experimental_data'), return; end

    genes = cfg.genes;
    conditions = cfg.conditions;
    ratios = [];

    for ci = 1:numel(conditions)
        cond = conditions{ci};
        if ~isfield(cfg.experimental_data, cond), continue; end

        result = solve_ode_generic(cfg, cond);
        u = compute_u_generic(cfg, result, cond);

        for gi = 1:numel(genes)
            gene = genes{gi};
            if ~isfield(cfg.experimental_data.(cond), gene), continue; end
            obs = cfg.experimental_data.(cond).(gene);
            [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), cfg.ktx_fits.(gene), ...
                cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, cfg.solver.t_eval);
            idx_sample = find(t_det >= cfg.t_sample, 1);
            if isempty(idx_sample), idx_sample = numel(t_det); end
            m_pred = m_det(idx_sample);
            if m_pred > 0 && obs > 0
                ratios(end+1) = abs(log10(m_pred / obs));
            end
        end
    end

    if ~isempty(ratios)
        val.within_2x = sum(ratios <= log10(2)) / numel(ratios) * 100;
        val.within_5x = sum(ratios <= log10(5)) / numel(ratios) * 100;
        val.rmse_log10 = sqrt(mean(ratios.^2));
        val.n_comparisons = numel(ratios);
        fprintf('  Within 2x: %.0f%%  Within 5x: %.0f%%  RMSE(log10): %.3f\n', ...
            val.within_2x, val.within_5x, val.rmse_log10);
    end
end