function [rankings, evaluations] = hypothesis_evaluator_generic(cfg, output_dir)
%HYPOTHESIS_EVALUATOR_GENERIC Universal 5-component hypothesis scoring.
%   [rankings, evaluations] = HYPOTHESIS_EVALUATOR_GENERIC(cfg, output_dir)
%
%   Scoring components:
%     C1 (0.25): Phenotype R2 — bridge prediction quality
%     C2 (0.30): Directional consistency — predicted up/down vs diagnostic_up
%     C3 (0.15): Magnitude accuracy — pred vs obs expression (if available)
%     C4 (0.20): Genomic support — fraction of diagnostic genes present
%     C5 (0.10): Parsimony — fewer diagnostic genes = more parsimonious

out_hyp = fullfile(output_dir, 'hypothesis_ranking');
if ~exist(out_hyp, 'dir'), mkdir(out_hyp); end

%% WEIGHTS
w = struct('phenotype', 0.25, 'directional', 0.30, 'magnitude', 0.15, ...
    'genomic', 0.20, 'parsimony', 0.10);

conditions = cfg.conditions;
genes      = cfg.genes;
n_genes    = numel(genes);

%% COMPUTE MODEL PREDICTIONS (per condition)
model_preds = struct();
for ci = 1:numel(conditions)
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    u = compute_u_generic(cfg, result, cond);
    for gi = 1:n_genes
        gene = genes{gi};
        ktx = cfg.ktx_fits.(gene);
        [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
            cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, cfg.solver.t_eval);
        idx_sample = find(t_det >= cfg.t_sample, 1);
        if isempty(idx_sample), idx_sample = numel(t_det); end
        model_preds.(cond).(gene) = m_det(idx_sample);
    end
end

%% COMPUTE BRIDGE R2 per condition (C1)
bridge_R2 = struct();
for ci = 1:numel(conditions)
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    phi = extract_phi_generic(cfg, result, cond);
    t = result.t(:)';
    D_ref = phenotype_reference_generic(cfg, cond, struct(), t);
    [~, ~, ~, R2, ~, ~] = fit_bridge(phi, t, D_ref);
    bridge_R2.(cond) = R2;
end

%% EVALUATE HYPOTHESES
evaluations = {};
idx = 0;

for ci = 1:numel(conditions)
    cond = conditions{ci};
    if ~isfield(cfg.hypotheses, cond), continue; end
    cond_hyps = cfg.hypotheses.(cond);
    for hi = 1:numel(cond_hyps)
        hyp = cond_hyps{hi};
        idx = idx + 1;
        ev = evaluate_single_hypothesis(hyp, cond, model_preds, bridge_R2, w, cfg);
        evaluations{idx} = ev; %#ok
    end
end

%% RANK PER CONDITION
rankings = struct();
for ci = 1:numel(conditions)
    cond = conditions{ci};
    cond_evs = {};
    for i = 1:numel(evaluations)
        if strcmp(evaluations{i}.condition, cond)
            cond_evs{end+1} = evaluations{i}; %#ok
        end
    end
    if isempty(cond_evs), continue; end
    scores = cellfun(@(e) e.score, cond_evs);
    [~, order] = sort(scores, 'descend');
    cond_evs = cond_evs(order);
    rankings.(cond).best = cond_evs{1};
    rankings.(cond).all  = cond_evs;
end

%% EXPORT
export_hypothesis_results(rankings, evaluations, conditions, out_hyp);

fprintf('=== Hypothesis evaluation complete ===\n');
for ci = 1:numel(conditions)
    cond = conditions{ci};
    if isfield(rankings, cond)
        fprintf('  %s winner: %s (score=%.4f)\n', cond, ...
            rankings.(cond).best.id, rankings.(cond).best.score);
    end
end
end

%% ═══ EVALUATE SINGLE HYPOTHESIS ═══
function ev = evaluate_single_hypothesis(hyp, cond, model_preds, bridge_R2, w, cfg)
    ev.id        = hyp.id;
    ev.label     = hyp.label;
    ev.condition = cond;

    %% C1: Phenotype R2 (bridge quality)
    ev.c1 = max(0, bridge_R2.(cond));

    %% C2: Directional consistency
    % Check if genes in diagnostic_up are actually predicted UP vs baseline
    diagnostic_up = hyp.diagnostic_up;
    % Use first condition as baseline for fold-change
    all_conds = fieldnames(model_preds);
    baseline_cond = all_conds{1};

    n_match = 0; n_check = 0;
    for gi = 1:numel(cfg.genes)
        gene = cfg.genes{gi};
        pred_val = model_preds.(cond).(gene);
        base_val = model_preds.(baseline_cond).(gene);
        is_diagnostic = any(strcmp(gene, diagnostic_up));

        if strcmp(cond, baseline_cond)
            % For the baseline condition, check if diagnostic genes have higher signal
            if is_diagnostic
                n_match = n_match + 1;
            end
            n_check = n_check + 1;
        else
            % For non-baseline: diagnostic genes should be upregulated
            fc = pred_val / max(base_val, 1e-12);
            if is_diagnostic
                if fc > 1.0, n_match = n_match + 1; end
            else
                if fc <= 1.2, n_match = n_match + 1; end  % non-diagnostic should be stable/down
            end
            n_check = n_check + 1;
        end
    end

    % Also integrate experimental data if available
    if isfield(cfg, 'experimental_data') && isfield(cfg.experimental_data, cond) && ...
       isfield(cfg.experimental_data, baseline_cond)
        for gi = 1:numel(cfg.genes)
            gene = cfg.genes{gi};
            if isfield(cfg.experimental_data.(cond), gene) && ...
               isfield(cfg.experimental_data.(baseline_cond), gene)
                obs_val = cfg.experimental_data.(cond).(gene);
                obs_base = cfg.experimental_data.(baseline_cond).(gene);
                obs_fc = obs_val / max(obs_base, 1e-12);

                pred_val = model_preds.(cond).(gene);
                base_val = model_preds.(baseline_cond).(gene);
                pred_fc = pred_val / max(base_val, 1e-12);

                % Direction match?
                if (obs_fc > 1 && pred_fc > 1) || (obs_fc < 1 && pred_fc < 1) || ...
                   (abs(obs_fc - 1) < 0.2 && abs(pred_fc - 1) < 0.2)
                    n_match = n_match + 1;
                end
                n_check = n_check + 1;
            end
        end
    end

    ev.c2 = n_match / max(n_check, 1);

    %% C3: Magnitude accuracy
    log_errors = [];
    if isfield(cfg, 'experimental_data') && isfield(cfg.experimental_data, cond)
        for gi = 1:numel(cfg.genes)
            gene = cfg.genes{gi};
            if isfield(cfg.experimental_data.(cond), gene)
                obs = cfg.experimental_data.(cond).(gene);
                pred = model_preds.(cond).(gene);
                if pred > 1e-12 && obs > 1e-12
                    log_errors(end+1) = abs(log10(pred / obs)); %#ok
                end
            end
        end
    end
    if ~isempty(log_errors)
        ev.c3 = max(0, 1 - mean(log_errors) / 3);
    else
        ev.c3 = 0.5;  % no data → neutral
    end

    %% C4: Genomic support
    % Fraction of diagnostic genes that are actually in the model genes list
    n_in_model = 0;
    for di = 1:numel(diagnostic_up)
        if any(strcmp(diagnostic_up{di}, cfg.genes))
            n_in_model = n_in_model + 1;
        end
    end
    ev.c4 = n_in_model / max(numel(diagnostic_up), 1);

    %% C5: Parsimony (fewer modules = more parsimonious)
    n_diag = numel(diagnostic_up);
    ev.c5 = 1 - (n_diag - 1) / max(numel(cfg.genes) - 1, 1);  % normalized

    %% PENALTY
    ev.penalty = 1.0;
    if ev.c2 < 0.4, ev.penalty = 0.5; end

    %% TOTAL SCORE
    ev.score = ev.penalty * (w.phenotype * ev.c1 + w.directional * ev.c2 + ...
        w.magnitude * ev.c3 + w.genomic * ev.c4 + w.parsimony * ev.c5);

    %% COMPARE WITH FROZEN SCORE
    if isfield(hyp, 'score_total')
        ev.frozen_score = hyp.score_total;
        ev.delta_frozen = ev.score - hyp.score_total;
    end

    ev.status = 'EVALUATED';
end

%% ═══ EXPORT ═══
function export_hypothesis_results(rankings, evaluations, conditions, out_dir)
    % per_condition_summary.tsv
    fid = fopen(fullfile(out_dir, 'per_condition_summary.tsv'), 'w');
    fprintf(fid, 'condition\tbest_hypothesis\tscore\tn_evaluated\n');
    for ci = 1:numel(conditions)
        cond = conditions{ci};
        if ~isfield(rankings, cond), continue; end
        r = rankings.(cond);
        fprintf(fid, '%s\t%s\t%.4f\t%d\n', cond, r.best.id, r.best.score, numel(r.all));
    end
    fclose(fid);

    % hypothesis_scores.tsv
    fid = fopen(fullfile(out_dir, 'hypothesis_scores.tsv'), 'w');
    fprintf(fid, 'hypothesis\tcondition\tscore\tc1_pheno\tc2_dir\tc3_mag\tc4_gen\tc5_pars\tpenalty\tstatus');
    if ~isempty(evaluations) && isfield(evaluations{1}, 'frozen_score')
        fprintf(fid, '\tfrozen_score\tdelta');
    end
    fprintf(fid, '\n');
    for i = 1:numel(evaluations)
        e = evaluations{i};
        fprintf(fid, '%s\t%s\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.2f\t%s', ...
            e.id, e.condition, e.score, e.c1, e.c2, e.c3, e.c4, e.c5, e.penalty, e.status);
        if isfield(e, 'frozen_score')
            fprintf(fid, '\t%.4f\t%.4f', e.frozen_score, e.delta_frozen);
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end