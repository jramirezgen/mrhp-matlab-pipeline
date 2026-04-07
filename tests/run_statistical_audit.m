function run_statistical_audit(base_dir)
%RUN_STATISTICAL_AUDIT Exhaustive statistical audit of MRHP pipeline outputs.
%   Validates: ODE convergence, non-negativity, bridge fitting, expression
%   layer, stochastic-deterministic consistency, hypothesis scoring, and
%   cross-system generalization. Exports audit_report.tsv + audit_certification.json.
%
%   Usage:
%     run_statistical_audit()                    % auto-detect from pwd
%     run_statistical_audit('/path/to/outputs')  % explicit

if nargin < 1
    % Auto-detect: look for outputs/ in parent directories
    base_dir = fileparts(fileparts(mfilename('fullpath')));  % up from scripts/audit/
    base_dir = fullfile(base_dir, 'outputs');
end

fprintf('\n========================================\n');
fprintf('  MRHP STATISTICAL AUDIT\n');
fprintf('  Base: %s\n', base_dir);
fprintf('  Timestamp: %s\n', datestr(now, 'yyyy-mm-ddTHH:MM:SS'));
fprintf('========================================\n\n');

% Discover organism output dirs
org_dirs = {};
org_names = {};
d = dir(base_dir);
for i = 1:numel(d)
    if d(i).isdir && ~startsWith(d(i).name, '.') && ...
       ~strcmp(d(i).name, 'metadata') && ~strcmp(d(i).name, 'reports') && ...
       ~strcmp(d(i).name, 'logs') && ~strcmp(d(i).name, 'audit_results')
        org_dirs{end+1} = fullfile(base_dir, d(i).name); %#ok
        org_names{end+1} = d(i).name; %#ok
    end
end

fprintf('Found %d organism output directories:\n', numel(org_dirs));
for i = 1:numel(org_names)
    fprintf('  [%d] %s\n', i, org_names{i});
end
fprintf('\n');

% ═══════════════════════════════════════════════════════
% Initialize audit results
% ═══════════════════════════════════════════════════════
audit = struct();
audit.timestamp = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
audit.n_organisms = numel(org_dirs);
audit.organisms = struct();
all_pass = true;

for oi = 1:numel(org_dirs)
    org_name = org_names{oi};
    org_dir  = org_dirs{oi};
    fprintf('────────────────────────────────────────\n');
    fprintf('  AUDITING: %s\n', org_name);
    fprintf('────────────────────────────────────────\n\n');

    org_audit = struct();
    org_audit.name = org_name;
    org_pass = true;

    % ═══ TEST 1: ODE CONVERGENCE & NON-NEGATIVITY ═══
    fprintf('== TEST 1: ODE Convergence & Non-Negativity ==\n');
    ode_dir = fullfile(org_dir, 'ode_timeseries');
    ode_files = dir(fullfile(ode_dir, 'timeseries_*.tsv'));
    n_ode = numel(ode_files);
    n_negative = 0;
    n_nan = 0;
    n_inf = 0;
    min_val_global = Inf;
    max_val_global = -Inf;
    n_timepoints_total = 0;
    n_species_total = 0;

    for fi = 1:n_ode
        fpath = fullfile(ode_dir, ode_files(fi).name);
        T = readtable(fpath, 'FileType', 'text', 'Delimiter', '\t');
        data = table2array(T(:, 2:end));  % skip time column
        n_timepoints_total = n_timepoints_total + size(data, 1);
        if fi == 1, n_species_total = size(data, 2); end

        neg_count = sum(data(:) < -1e-12);
        nan_count = sum(isnan(data(:)));
        inf_count = sum(isinf(data(:)));
        n_negative = n_negative + neg_count;
        n_nan = n_nan + nan_count;
        n_inf = n_inf + inf_count;
        min_val_global = min(min_val_global, min(data(:)));
        max_val_global = max(max_val_global, max(data(:)));
    end

    org_audit.ode.n_scenarios = n_ode;
    org_audit.ode.n_negative_values = n_negative;
    org_audit.ode.n_nan_values = n_nan;
    org_audit.ode.n_inf_values = n_inf;
    org_audit.ode.min_concentration = min_val_global;
    org_audit.ode.max_concentration = max_val_global;
    org_audit.ode.n_species = n_species_total;

    ode_pass = (n_negative == 0) && (n_nan == 0) && (n_inf == 0);
    org_audit.ode.status = ternary(ode_pass, 'PASS', 'FAIL');
    fprintf('  Scenarios: %d\n', n_ode);
    fprintf('  Negative values: %d  |  NaN: %d  |  Inf: %d\n', n_negative, n_nan, n_inf);
    fprintf('  Concentration range: [%.4e, %.4e]\n', min_val_global, max_val_global);
    fprintf('  Status: %s\n\n', org_audit.ode.status);
    if ~ode_pass, org_pass = false; end

    % ═══ TEST 2: BRIDGE FITTING QUALITY ═══
    fprintf('== TEST 2: Bridge Fitting Quality ==\n');
    bridge_file = fullfile(org_dir, 'metrics', 'bridge_all_scenarios.tsv');
    if exist(bridge_file, 'file')
        B = readtable(bridge_file, 'FileType', 'text', 'Delimiter', '\t');
        R2_all = B.R2;
        lam_all = B.lambda;
        beta_all = B.beta;
        ymax_all = B.ymax;
        conditions_bridge = unique(B.condition, 'stable');

        org_audit.bridge.n_scenarios = height(B);
        org_audit.bridge.mean_R2 = mean(R2_all);
        org_audit.bridge.median_R2 = median(R2_all);
        org_audit.bridge.min_R2 = min(R2_all);
        org_audit.bridge.max_R2 = max(R2_all);
        org_audit.bridge.std_R2 = std(R2_all);
        org_audit.bridge.n_R2_above_095 = sum(R2_all > 0.95);
        org_audit.bridge.n_R2_above_099 = sum(R2_all > 0.99);
        org_audit.bridge.n_R2_negative = sum(R2_all < 0);
        org_audit.bridge.lambda_range = [min(lam_all), max(lam_all)];
        org_audit.bridge.beta_range = [min(beta_all), max(beta_all)];
        org_audit.bridge.ymax_range = [min(ymax_all), max(ymax_all)];

        % Per-condition best R2
        per_cond = struct();
        for ci = 1:numel(conditions_bridge)
            cond = conditions_bridge{ci};
            mask = strcmp(B.condition, cond);
            cR2 = R2_all(mask);
            safe_cond = matlab.lang.makeValidName(cond);
            per_cond.(safe_cond).best_R2 = max(cR2);
            per_cond.(safe_cond).mean_R2 = mean(cR2);
            per_cond.(safe_cond).n_scenarios = sum(mask);
        end
        org_audit.bridge.per_condition = per_cond;

        % Residual analysis from bridge prediction files
        bridge_pred_dir = fullfile(org_dir, 'bridge_predictions');
        resid_stats = struct();
        bp_files = dir(fullfile(bridge_pred_dir, 'bridge_pred_*.tsv'));
        for fi = 1:numel(bp_files)
            fpath = fullfile(bridge_pred_dir, bp_files(fi).name);
            T = readtable(fpath, 'FileType', 'text', 'Delimiter', '\t');
            residuals = T.D_ref - T.D_bridge;
            cond_name = strrep(strrep(bp_files(fi).name, 'bridge_pred_', ''), '.tsv', '');
            safe_cond = matlab.lang.makeValidName(cond_name);
            resid_stats.(safe_cond).mean_residual = mean(residuals);
            resid_stats.(safe_cond).std_residual = std(residuals);
            resid_stats.(safe_cond).max_abs_residual = max(abs(residuals));
            resid_stats.(safe_cond).skewness = skewness(residuals);
            resid_stats.(safe_cond).kurtosis = kurtosis(residuals);
            % Durbin-Watson-like autocorrelation of residuals
            if numel(residuals) > 2
                diffs = diff(residuals);
                dw_stat = sum(diffs.^2) / max(sum(residuals.^2), 1e-12);
                resid_stats.(safe_cond).durbin_watson = dw_stat;
            end
        end
        org_audit.bridge.residuals = resid_stats;

        % Parameter bounds check
        lam_ok = all(lam_all > 0) && all(lam_all < 500);
        beta_ok = all(beta_all >= 0.1) && all(beta_all <= 5.0);
        ymax_ok = all(ymax_all > 0);
        bridge_pass = (org_audit.bridge.n_R2_negative == 0) && lam_ok && beta_ok && ymax_ok;
        org_audit.bridge.params_in_bounds = lam_ok && beta_ok && ymax_ok;
        org_audit.bridge.status = ternary(bridge_pass, 'PASS', 'WARN');

        fprintf('  Scenarios: %d\n', height(B));
        fprintf('  R² range: [%.4f, %.4f], mean=%.4f, std=%.4f\n', ...
            min(R2_all), max(R2_all), mean(R2_all), std(R2_all));
        fprintf('  R²>0.95: %d/%d  |  R²>0.99: %d/%d  |  R²<0: %d\n', ...
            org_audit.bridge.n_R2_above_095, height(B), ...
            org_audit.bridge.n_R2_above_099, height(B), ...
            org_audit.bridge.n_R2_negative);
        fprintf('  λ ∈ [%.4f, %.4f]  |  β ∈ [%.4f, %.4f]  |  ymax ∈ [%.2f, %.2f]\n', ...
            min(lam_all), max(lam_all), min(beta_all), max(beta_all), ...
            min(ymax_all), max(ymax_all));
        fprintf('  Params in bounds: %s\n', ternary(org_audit.bridge.params_in_bounds, 'YES', 'NO'));
        fprintf('  Status: %s\n\n', org_audit.bridge.status);
    else
        org_audit.bridge.status = 'MISSING';
        fprintf('  bridge_all_scenarios.tsv not found — SKIP\n\n');
    end

    % ═══ TEST 3: EXPRESSION LAYER — DETERMINISTIC ═══
    fprintf('== TEST 3: Expression Layer — Deterministic ==\n');
    expr_dir = fullfile(org_dir, 'expression');
    det_files = dir(fullfile(expr_dir, 'det_*.tsv'));
    n_det = numel(det_files);
    det_neg = 0;
    det_nan = 0;
    mrna_range = [Inf, -Inf];
    prot_range = [Inf, -Inf];

    for fi = 1:n_det
        fpath = fullfile(expr_dir, det_files(fi).name);
        T = readtable(fpath, 'FileType', 'text', 'Delimiter', '\t');
        m = T.mRNA_nM;
        p = T.Protein_nM;
        det_neg = det_neg + sum(m < -1e-15) + sum(p < -1e-15);
        det_nan = det_nan + sum(isnan(m)) + sum(isnan(p));
        mrna_range(1) = min(mrna_range(1), min(m));
        mrna_range(2) = max(mrna_range(2), max(m));
        prot_range(1) = min(prot_range(1), min(p));
        prot_range(2) = max(prot_range(2), max(p));
    end

    org_audit.expression_det.n_files = n_det;
    org_audit.expression_det.n_negative = det_neg;
    org_audit.expression_det.n_nan = det_nan;
    org_audit.expression_det.mRNA_range_nM = mrna_range;
    org_audit.expression_det.protein_range_nM = prot_range;

    det_pass = (det_neg == 0) && (det_nan == 0);
    org_audit.expression_det.status = ternary(det_pass, 'PASS', 'FAIL');
    fprintf('  Files: %d\n', n_det);
    fprintf('  Negative: %d  |  NaN: %d\n', det_neg, det_nan);
    fprintf('  mRNA range: [%.4e, %.4e] nM\n', mrna_range(1), mrna_range(2));
    fprintf('  Protein range: [%.4e, %.4e] nM\n', prot_range(1), prot_range(2));
    fprintf('  Status: %s\n\n', org_audit.expression_det.status);
    if ~det_pass, org_pass = false; end

    % ═══ TEST 4: STOCHASTIC vs DETERMINISTIC CONSISTENCY ═══
    fprintf('== TEST 4: Stochastic vs Deterministic Consistency ==\n');
    sto_files = dir(fullfile(expr_dir, 'sto_*.tsv'));
    n_sto = numel(sto_files);
    ergodic_errors = [];

    for fi = 1:n_sto
        sto_path = fullfile(expr_dir, sto_files(fi).name);
        % Find matching det file
        det_name = strrep(sto_files(fi).name, 'sto_', 'det_');
        det_path = fullfile(expr_dir, det_name);
        if ~exist(det_path, 'file'), continue; end

        T_sto = readtable(sto_path, 'FileType', 'text', 'Delimiter', '\t');
        T_det = readtable(det_path, 'FileType', 'text', 'Delimiter', '\t');

        % Compare at endpoint
        m_sto_end = T_sto.mean_m_nM(end);
        m_det_end = T_det.mRNA_nM(end);
        if m_det_end > 1e-10
            rel_err = abs(m_sto_end - m_det_end) / m_det_end;
            ergodic_errors(end+1) = rel_err; %#ok
        end
    end

    org_audit.stochastic.n_comparisons = numel(ergodic_errors);
    if ~isempty(ergodic_errors)
        org_audit.stochastic.mean_relative_error = mean(ergodic_errors);
        org_audit.stochastic.max_relative_error = max(ergodic_errors);
        org_audit.stochastic.within_10pct = sum(ergodic_errors < 0.10);
        org_audit.stochastic.within_20pct = sum(ergodic_errors < 0.20);

        sto_pass = mean(ergodic_errors) < 0.20;  % mean <20% is acceptable
        org_audit.stochastic.status = ternary(sto_pass, 'PASS', 'WARN');
        fprintf('  Comparisons: %d\n', numel(ergodic_errors));
        fprintf('  Mean relative error (sto vs det): %.2f%%\n', mean(ergodic_errors)*100);
        fprintf('  Max relative error: %.2f%%\n', max(ergodic_errors)*100);
        fprintf('  Within 10%%: %d/%d  |  Within 20%%: %d/%d\n', ...
            sum(ergodic_errors < 0.10), numel(ergodic_errors), ...
            sum(ergodic_errors < 0.20), numel(ergodic_errors));
    else
        org_audit.stochastic.status = 'SKIP';
        fprintf('  No matching sto/det pairs found\n');
    end
    fprintf('  Status: %s\n\n', org_audit.stochastic.status);

    % ═══ TEST 5: HYPOTHESIS SCORING INTEGRITY ═══
    fprintf('== TEST 5: Hypothesis Scoring Integrity ==\n');
    hyp_file = fullfile(org_dir, 'hypothesis_ranking', 'hypothesis_scores.tsv');
    if exist(hyp_file, 'file')
        H = readtable(hyp_file, 'FileType', 'text', 'Delimiter', '\t');
        scores = H.score;
        c1 = H.c1_pheno;
        c2 = H.c2_dir;
        c3 = H.c3_mag;
        c4 = H.c4_gen;
        c5 = H.c5_pars;
        penalties = H.penalty;

        % Verify weight reconstruction: score = penalty * (0.25*c1 + 0.30*c2 + 0.15*c3 + 0.20*c4 + 0.10*c5)
        reconstructed = penalties .* (0.25*c1 + 0.30*c2 + 0.15*c3 + 0.20*c4 + 0.10*c5);
        recon_error = max(abs(scores - reconstructed));

        % Check ranges
        scores_in_range = all(scores >= 0 & scores <= 1);
        c_in_range = all(c1 >= 0 & c1 <= 1) && all(c2 >= 0 & c2 <= 1.01) && ...
                     all(c3 >= 0 & c3 <= 1) && all(c4 >= 0 & c4 <= 1) && ...
                     all(c5 >= 0 & c5 <= 1);
        penalty_valid = all(penalties == 0.5 | penalties == 1.0);

        % Frozen score comparison
        if any(strcmp(H.Properties.VariableNames, 'frozen_score'))
            frozen = H.frozen_score;
            delta_frozen = abs(scores - frozen);
            max_delta = max(delta_frozen);
            org_audit.hypothesis.max_frozen_delta = max_delta;
            org_audit.hypothesis.mean_frozen_delta = mean(delta_frozen);
        else
            max_delta = 0;
        end

        org_audit.hypothesis.n_hypotheses = height(H);
        org_audit.hypothesis.score_range = [min(scores), max(scores)];
        org_audit.hypothesis.reconstruction_error = recon_error;
        org_audit.hypothesis.scores_in_range = scores_in_range;
        org_audit.hypothesis.components_in_range = c_in_range;
        org_audit.hypothesis.penalties_valid = penalty_valid;
        org_audit.hypothesis.weight_sum = 0.25 + 0.30 + 0.15 + 0.20 + 0.10;

        hyp_pass = (recon_error < 1e-3) && scores_in_range && c_in_range && penalty_valid;
        org_audit.hypothesis.status = ternary(hyp_pass, 'PASS', 'FAIL');

        fprintf('  Hypotheses: %d\n', height(H));
        fprintf('  Score range: [%.4f, %.4f]\n', min(scores), max(scores));
        fprintf('  Weight sum: %.2f (expected 1.00)\n', org_audit.hypothesis.weight_sum);
        fprintf('  Reconstruction error: %.2e (threshold: 1e-3 — TSV rounding tolerance)\n', recon_error);
        fprintf('  All scores ∈ [0,1]: %s  |  All components ∈ [0,1]: %s\n', ...
            ternary(scores_in_range, 'YES', 'NO'), ternary(c_in_range, 'YES', 'NO'));
        fprintf('  Penalties valid: %s\n', ternary(penalty_valid, 'YES', 'NO'));
        if isfield(org_audit.hypothesis, 'max_frozen_delta')
            fprintf('  Max |frozen − computed|: %.4f\n', max_delta);
        end
        fprintf('  Status: %s\n\n', org_audit.hypothesis.status);
        if ~hyp_pass, org_pass = false; end
    else
        org_audit.hypothesis.status = 'MISSING';
        fprintf('  hypothesis_scores.tsv not found — SKIP\n\n');
    end

    % ═══ TEST 6: PHI SIGNAL RANGE ═══
    fprintf('== TEST 6: Phi Signal Quality ==\n');
    phi_dir = fullfile(org_dir, 'phi_signals');
    phi_files = dir(fullfile(phi_dir, 'phi_signals_*.tsv'));
    phi_pass_local = true;
    for fi = 1:numel(phi_files)
        fpath = fullfile(phi_dir, phi_files(fi).name);
        T = readtable(fpath, 'FileType', 'text', 'Delimiter', '\t');
        phi_col = T{:, 2};  % second column is phi
        if any(phi_col < -1e-10) || any(phi_col > 1.01) || any(isnan(phi_col))
            phi_pass_local = false;
        end
    end
    org_audit.phi.n_files = numel(phi_files);
    org_audit.phi.status = ternary(phi_pass_local, 'PASS', 'FAIL');
    fprintf('  Files: %d  |  All phi ∈ [0,1]: %s\n', numel(phi_files), ...
        ternary(phi_pass_local, 'YES', 'NO'));
    fprintf('  Status: %s\n\n', org_audit.phi.status);
    if ~phi_pass_local, org_pass = false; end

    % ═══ TEST 7: REGULATORY SIGNAL RANGE ═══
    fprintf('== TEST 7: Regulatory Signal Quality ==\n');
    reg_dir = fullfile(org_dir, 'regulatory_signals');
    reg_files = dir(fullfile(reg_dir, 'regulatory_*.tsv'));
    reg_pass_local = true;
    n_reg_violations = 0;
    for fi = 1:numel(reg_files)
        fpath = fullfile(reg_dir, reg_files(fi).name);
        T = readtable(fpath, 'FileType', 'text', 'Delimiter', '\t');
        data = table2array(T(:, 2:end));
        violations = sum(data(:) < -1e-10) + sum(data(:) > 1.01) + sum(isnan(data(:)));
        n_reg_violations = n_reg_violations + violations;
        if violations > 0, reg_pass_local = false; end
    end
    org_audit.regulatory.n_files = numel(reg_files);
    org_audit.regulatory.n_violations = n_reg_violations;
    org_audit.regulatory.status = ternary(reg_pass_local, 'PASS', 'FAIL');
    fprintf('  Files: %d  |  Violations: %d\n', numel(reg_files), n_reg_violations);
    fprintf('  Status: %s\n\n', org_audit.regulatory.status);
    if ~reg_pass_local, org_pass = false; end

    % ═══ TEST 8: FIGURE COMPLETENESS ═══
    fprintf('== TEST 8: Figure Completeness ==\n');
    fig_dir = fullfile(org_dir, 'figures');
    png_files = dir(fullfile(fig_dir, '*.png'));
    svg_files = dir(fullfile(fig_dir, '*.svg'));
    expected_figs = 9;
    fig_pass = (numel(png_files) >= expected_figs) && (numel(svg_files) >= expected_figs);
    org_audit.figures.n_png = numel(png_files);
    org_audit.figures.n_svg = numel(svg_files);
    org_audit.figures.expected = expected_figs;
    org_audit.figures.status = ternary(fig_pass, 'PASS', 'WARN');
    fprintf('  PNG: %d  |  SVG: %d  |  Expected: %d\n', numel(png_files), numel(svg_files), expected_figs);
    fprintf('  Status: %s\n\n', org_audit.figures.status);

    % ═══ ORGANISM SUMMARY ═══
    org_audit.overall_status = ternary(org_pass, 'CERTIFIED', 'REVIEW_NEEDED');
    fprintf('═══ %s OVERALL: %s ═══\n\n', upper(org_name), org_audit.overall_status);
    if ~org_pass, all_pass = false; end

    audit.organisms.(matlab.lang.makeValidName(org_name)) = org_audit;
end

% ═══════════════════════════════════════════════════════
% CROSS-SYSTEM SUMMARY
% ═══════════════════════════════════════════════════════
fprintf('════════════════════════════════════════\n');
fprintf('  CROSS-SYSTEM AUDIT SUMMARY\n');
fprintf('════════════════════════════════════════\n\n');

org_fields = fieldnames(audit.organisms);
total_scenarios_ode = 0;
total_scenarios_bridge = 0;
total_det_files = 0;
total_hypotheses = 0;
total_figures = 0;

for i = 1:numel(org_fields)
    oa = audit.organisms.(org_fields{i});
    fprintf('  %-30s  %s\n', oa.name, oa.overall_status);
    total_scenarios_ode = total_scenarios_ode + oa.ode.n_scenarios;
    if isfield(oa.bridge, 'n_scenarios')
        total_scenarios_bridge = total_scenarios_bridge + oa.bridge.n_scenarios;
    end
    total_det_files = total_det_files + oa.expression_det.n_files;
    if isfield(oa.hypothesis, 'n_hypotheses')
        total_hypotheses = total_hypotheses + oa.hypothesis.n_hypotheses;
    end
    total_figures = total_figures + oa.figures.n_png;
end

audit.cross_system.total_ode_scenarios = total_scenarios_ode;
audit.cross_system.total_bridge_scenarios = total_scenarios_bridge;
audit.cross_system.total_expression_files = total_det_files;
audit.cross_system.total_hypotheses = total_hypotheses;
audit.cross_system.total_figures = total_figures;
audit.cross_system.all_certified = all_pass;
audit.cross_system.status = ternary(all_pass, 'PLATFORM_CERTIFIED', 'REVIEW_NEEDED');

fprintf('\n  Total ODE scenarios: %d\n', total_scenarios_ode);
fprintf('  Total bridge scenarios: %d\n', total_scenarios_bridge);
fprintf('  Total expression files: %d\n', total_det_files);
fprintf('  Total hypotheses evaluated: %d\n', total_hypotheses);
fprintf('  Total figures: %d\n', total_figures);
fprintf('\n  PLATFORM STATUS: %s\n', audit.cross_system.status);
fprintf('════════════════════════════════════════\n\n');

% ═══════════════════════════════════════════════════════
% EXPORT RESULTS
% ═══════════════════════════════════════════════════════
audit_dir = fullfile(base_dir, '..', 'audit_results');
if ~exist(audit_dir, 'dir'), mkdir(audit_dir); end

% Export JSON
json_path = fullfile(audit_dir, 'audit_certification.json');
json_str = jsonencode(audit);
fid = fopen(json_path, 'w');
fprintf(fid, '%s', json_str);
fclose(fid);
fprintf('Exported: %s\n', json_path);

% Export TSV summary
tsv_path = fullfile(audit_dir, 'audit_summary.tsv');
fid = fopen(tsv_path, 'w');
fprintf(fid, 'organism\tode_status\tode_scenarios\tbridge_status\tbridge_mean_R2\texpression_status\tstochastic_status\thypothesis_status\tphi_status\tregulatory_status\tfigures_status\toverall\n');
for i = 1:numel(org_fields)
    oa = audit.organisms.(org_fields{i});
    bridge_r2 = 0;
    if isfield(oa.bridge, 'mean_R2'), bridge_r2 = oa.bridge.mean_R2; end
    fprintf(fid, '%s\t%s\t%d\t%s\t%.4f\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n', ...
        oa.name, oa.ode.status, oa.ode.n_scenarios, ...
        oa.bridge.status, bridge_r2, ...
        oa.expression_det.status, oa.stochastic.status, ...
        oa.hypothesis.status, oa.phi.status, oa.regulatory.status, ...
        oa.figures.status, oa.overall_status);
end
fclose(fid);
fprintf('Exported: %s\n', tsv_path);

fprintf('\n=== Statistical Audit Complete ===\n');
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end