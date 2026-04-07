function figure_generation_generic(cfg, output_dir)
%FIGURE_GENERATION_GENERIC Generate publication figures for any organism.
%   Generates ~10 core figures adaptively based on config dimensions.

if nargin < 2
    output_dir = cfg.output_dir;
end

fig_dir = fullfile(output_dir, 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

conditions = cfg.conditions;
genes      = cfg.genes;
n_cond     = numel(conditions);
n_genes    = numel(genes);
t_eval     = cfg.solver.t_eval;

is_coupled = isfield(cfg, 'coupled_expression') && cfg.coupled_expression;

fprintf('=== Figure Generation ===\n\n');

%% FIG 1: ODE TIMESERIES
fprintf('  FIG 1: ODE timeseries ...\n');
n_cols = min(n_cond, 4);
n_rows = ceil(n_cond / n_cols);
fig1 = figure('Position',[100 100 400*n_cols 400*n_rows],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    t = result.t;
    model = cfg.models.(cond);
    ax = subplot(n_rows, n_cols, ci);
    hold(ax, 'on');
    % Plot target species
    if isfield(model, 'target_idx')
        plot(ax, t, result.y(model.target_idx,:), 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
    end
    % Plot biomass
    if isfield(model, 'biomass_idx')
        plot(ax, t, result.y(model.biomass_idx,:), 'Color', '#2196F3', 'LineWidth', 1.5);
    end
    xlabel(ax, 'Time (h)'); ylabel(ax, 'Concentration (mM)');
    title(ax, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    xlim(ax, cfg.solver.t_span);
    if ci == 1
        leg_entries = {};
        if isfield(model, 'target_idx'), leg_entries{end+1} = 'Target'; end
        if isfield(model, 'biomass_idx'), leg_entries{end+1} = 'Biomass'; end
        legend(ax, leg_entries, 'Location', 'east');
    end
end
sgtitle(sprintf('%s — ODE Metabolic Models', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig1, fullfile(fig_dir, 'fig1_ode_timeseries.png'), 'Resolution', 300);
saveas(fig1, fullfile(fig_dir, 'fig1_ode_timeseries.svg'));
close(fig1);

%% FIG 2: REGULATORY SIGNALS
fprintf('  FIG 2: Regulatory signals ...\n');
n_gcols = min(n_genes, 3);
n_grows = ceil(n_genes / n_gcols);
fig2 = figure('Position',[100 100 400*n_gcols 300*n_grows],'Visible','off');
for gi = 1:n_genes
    gene = genes{gi};
    ax = subplot(n_grows, n_gcols, gi);
    hold(ax, 'on');
    for ci = 1:n_cond
        cond = conditions{ci};
        result = solve_ode_generic(cfg, cond);
        u = compute_u_generic(cfg, result, cond);
        if isfield(u, gene)
            plot(ax, u.t, u.(gene), 'Color', get_cond_color(cfg, cond), 'LineWidth', 1.5);
        end
    end
    title(ax, gene, 'FontWeight', 'bold');
    xlim(ax, cfg.solver.t_span); ylim(ax, [0 1]);
    xlabel(ax, 'Time (h)'); ylabel(ax, 'u(t)');
    if gi == 1, legend(ax, conditions, 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('%s — Regulatory Signals', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig2, fullfile(fig_dir, 'fig2_regulatory_signals.png'), 'Resolution', 300);
saveas(fig2, fullfile(fig_dir, 'fig2_regulatory_signals.svg'));
close(fig2);

%% FIG 3: EXPRESSION TIME COURSES
fprintf('  FIG 3: Expression time courses ...\n');
fig3 = figure('Position',[100 100 400*n_gcols 300*n_grows],'Visible','off');
for gi = 1:n_genes
    gene = genes{gi};
    ax = subplot(n_grows, n_gcols, gi);
    hold(ax, 'on');
    for ci = 1:n_cond
        cond = conditions{ci};
        result = solve_ode_generic(cfg, cond);
        if is_coupled
            m_det = result.y(18 + gi, :);
            plot(ax, result.t, m_det, 'Color', get_cond_color(cfg, cond), 'LineWidth', 1.5);
        else
            u = compute_u_generic(cfg, result, cond);
            if isfield(u, gene)
                ktx = cfg.ktx_fits.(gene);
                [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
                    cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
                plot(ax, t_det, m_det, 'Color', get_cond_color(cfg, cond), 'LineWidth', 1.5);
            end
        end
    end
    title(ax, gene, 'FontWeight', 'bold');
    xlim(ax, cfg.solver.t_span);
    xlabel(ax, 'Time (h)'); ylabel(ax, 'mRNA (nM)');
    if gi == 1, legend(ax, conditions, 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('%s — Expression Dynamics', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig3, fullfile(fig_dir, 'fig3_expression_timecourses.png'), 'Resolution', 300);
saveas(fig3, fullfile(fig_dir, 'fig3_expression_timecourses.svg'));
close(fig3);

%% FIG 4: BRIDGE FITTING
fprintf('  FIG 4: Bridge fitting ...\n');
fig4 = figure('Position',[100 100 400*n_cols 400*ceil(n_cond/n_cols)],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    ax = subplot(ceil(n_cond/n_cols), n_cols, ci);
    hold(ax, 'on');
    result = solve_ode_generic(cfg, cond);
    t = result.t(:)';
    phi = extract_phi_generic(cfg, result, cond);
    D_ref = phenotype_reference_generic(cfg, cond, struct(), t);
    [~, ~, ~, R2, ~, y_pred] = fit_bridge(phi, t, D_ref);
    plot(ax, t, D_ref, 'k--', 'LineWidth', 1.5);
    plot(ax, t, y_pred, 'Color', get_cond_color(cfg, cond), 'LineWidth', 2.5);
    xlabel(ax, 'Time (h)'); ylabel(ax, cfg.phenotype_label);
    title(ax, sprintf('%s  R^2=%.4f', get_cond_label(cfg, cond), R2), ...
          'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    xlim(ax, cfg.solver.t_span);
    if ci == 1, legend(ax, {'Reference','Bridge'}, 'Location', 'southeast'); end
end
sgtitle(sprintf('%s — Phenotypic Bridge', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig4, fullfile(fig_dir, 'fig4_bridge.png'), 'Resolution', 300);
saveas(fig4, fullfile(fig_dir, 'fig4_bridge.svg'));
close(fig4);

%% FIG 5: INTEGRATED MULTISCALE
fprintf('  FIG 5: Integrated multiscale ...\n');
n_rows_ms = min(4, 2 + min(n_genes, 2));
fig5 = figure('Position',[100 100 400*n_cols 250*n_rows_ms],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    t = result.t;
    u = compute_u_generic(cfg, result, cond);
    model = cfg.models.(cond);
    % Row 1: Metabolism
    ax1 = subplot(n_rows_ms, n_cond, ci);
    hold(ax1, 'on');
    if isfield(model, 'target_idx')
        plot(ax1, t, result.y(model.target_idx,:), 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
    end
    if isfield(model, 'biomass_idx')
        plot(ax1, t, result.y(model.biomass_idx,:), 'Color', '#2196F3', 'LineWidth', 1.5);
    end
    title(ax1, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    xlim(ax1, cfg.solver.t_span);
    if ci == 1, ylabel(ax1, 'Metabolism'); end
    % Row 2: Bridge
    ax2 = subplot(n_rows_ms, n_cond, n_cond + ci);
    phi = extract_phi_generic(cfg, result, cond);
    D_ref = phenotype_reference_generic(cfg, cond, struct(), t);
    [~,~,~,~,~, y_pred] = fit_bridge(phi, t, D_ref);
    plot(ax2, t, y_pred, 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
    xlim(ax2, cfg.solver.t_span);
    if ci == 1, ylabel(ax2, 'Bridge'); end
    % Rows 3+: Top genes expression (from coupled ODE or legacy)
    for gri = 1:min(n_genes, n_rows_ms-2)
        gene = genes{gri};
        ax_g = subplot(n_rows_ms, n_cond, (1+gri)*n_cond + ci);
        if is_coupled
            plot(ax_g, t, result.y(18 + gri, :), 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
        elseif isfield(u, gene)
            ktx = cfg.ktx_fits.(gene);
            [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
                cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
            plot(ax_g, t_det, m_det, 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
        end
        xlim(ax_g, cfg.solver.t_span);
        if ci == 1, ylabel(ax_g, sprintf('%s mRNA', gene)); end
        if gri == min(n_genes, n_rows_ms-2), xlabel(ax_g, 'Time (h)'); end
    end
end
sgtitle(sprintf('%s — Integrated Multiscale View', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig5, fullfile(fig_dir, 'fig5_integrated_multiscale.png'), 'Resolution', 300);
saveas(fig5, fullfile(fig_dir, 'fig5_integrated_multiscale.svg'));
close(fig5);

%% FIG 6: STOCHASTIC vs DETERMINISTIC
fprintf('  FIG 6: Stochastic vs Deterministic ...\n');
ref_cond = conditions{1};
result = solve_ode_generic(cfg, ref_cond);
u = compute_u_generic(cfg, result, ref_cond);
expr_params = struct('stochastic_tau', cfg.expression.stochastic_tau, ...
    'beta_m', cfg.expression.beta_m, 'beta_p', cfg.expression.beta_p, ...
    'ktl', cfg.expression.ktl, 'nm_to_mol', cfg.expression.nm_to_mol, ...
    't_span', cfg.solver.t_span, 't_eval', cfg.solver.t_eval);
fig6 = figure('Position',[100 100 400*n_gcols 300*n_grows],'Visible','off');
for gi = 1:n_genes
    gene = genes{gi};
    ax = subplot(n_grows, n_gcols, gi);
    hold(ax, 'on');
    if isfield(u, gene)
        ktx = cfg.ktx_fits.(gene);
        [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
            cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
        sto = tau_leaping_expression(u.t, u.(gene), ktx, cfg.expression.n_cells, 42, expr_params);
        fill(ax, [sto.t fliplr(sto.t)], [sto.mean_m+sto.std_m fliplr(sto.mean_m-sto.std_m)], ...
            get_cond_color(cfg, ref_cond), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        plot(ax, sto.t, sto.mean_m, 'Color', get_cond_color(cfg, ref_cond), 'LineWidth', 1.5);
        plot(ax, t_det, m_det, 'k--', 'LineWidth', 1);
    end
    title(ax, gene, 'FontWeight', 'bold');
    xlim(ax, cfg.solver.t_span);
    xlabel(ax, 'Time (h)'); ylabel(ax, 'mRNA (nM)');
    if gi == 1, legend(ax, {'±σ','Stochastic','Deterministic'}, 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('%s — Stochastic vs Deterministic (%s)', cfg.organism, ref_cond), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig6, fullfile(fig_dir, 'fig6_stochastic_vs_det.png'), 'Resolution', 300);
saveas(fig6, fullfile(fig_dir, 'fig6_stochastic_vs_det.svg'));
close(fig6);

%% FIG 7: HYPOTHESIS RANKING
fprintf('  FIG 7: Hypothesis ranking ...\n');
if isfield(cfg, 'hypotheses') && ~isempty(fieldnames(cfg.hypotheses))
    fig7 = figure('Position',[100 100 500*n_cols 400*ceil(n_cond/n_cols)],'Visible','off');
    for ci = 1:n_cond
        cond = conditions{ci};
        if ~isfield(cfg.hypotheses, cond), continue; end
        hyps = cfg.hypotheses.(cond);
        ax = subplot(ceil(n_cond/n_cols), n_cols, ci);
        names = {}; scores = [];
        for hi = 1:numel(hyps)
            h = hyps{hi};
            names{hi} = h.id;
            scores(hi) = h.score_total;
        end
        barh(ax, scores, 'FaceColor', get_cond_color(cfg, cond));
        set(ax, 'YTick', 1:numel(names), 'YTickLabel', strrep(names, '_', '\_'));
        xlabel(ax, 'Score');
        title(ax, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
        xlim(ax, [0 1]);
    end
    sgtitle(sprintf('%s — Hypothesis Ranking', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
    exportgraphics(fig7, fullfile(fig_dir, 'fig7_hypothesis_ranking.png'), 'Resolution', 300);
    saveas(fig7, fullfile(fig_dir, 'fig7_hypothesis_ranking.svg'));
    close(fig7);
else
    fprintf('    No hypotheses defined — skipping\n');
end

%% FIG 8: CONDITION COMPARISON
fprintf('  FIG 8: Condition comparison ...\n');
fig8 = figure('Position',[100 100 600 400],'Visible','off');
ax = axes; hold(ax, 'on');
m_at_sample = nan(n_cond, n_genes);
for ci = 1:n_cond
    cond = conditions{ci};
    result = solve_ode_generic(cfg, cond);
    u = compute_u_generic(cfg, result, cond);
    for gi = 1:n_genes
        gene = genes{gi};
        if is_coupled
            m_trace = result.y(18 + gi, :);
            idx_sample = find(result.t >= cfg.t_sample, 1);
            if isempty(idx_sample), idx_sample = numel(result.t); end
            m_at_sample(ci, gi) = m_trace(idx_sample);
        else
            if isfield(u, gene)
                ktx = cfg.ktx_fits.(gene);
                [t_det, m_det, ~] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
                    cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
                idx_sample = find(t_det >= cfg.t_sample, 1);
                if isempty(idx_sample), idx_sample = numel(t_det); end
                m_at_sample(ci, gi) = m_det(idx_sample);
            end
        end
    end
end
bar(ax, m_at_sample);
set(ax, 'XTick', 1:n_cond, 'XTickLabel', conditions);
legend(ax, genes, 'Location', 'northwest', 'FontSize', 7);
ylabel(ax, sprintf('mRNA @%.0fh (nM)', cfg.t_sample));
title(ax, sprintf('%s — Expression Comparison', cfg.organism), 'FontWeight', 'bold');
exportgraphics(fig8, fullfile(fig_dir, 'fig8_condition_comparison.png'), 'Resolution', 300);
saveas(fig8, fullfile(fig_dir, 'fig8_condition_comparison.svg'));
close(fig8);

%% FIG 9: ENZYME / PROTEIN DYNAMICS
fprintf('  FIG 9: Enzyme / Protein dynamics ...\n');
fig9 = figure('Position',[100 100 400*n_gcols 300*n_grows],'Visible','off');
for gi = 1:n_genes
    gene = genes{gi};
    ax = subplot(n_grows, n_gcols, gi);
    hold(ax, 'on');
    for ci = 1:n_cond
        cond = conditions{ci};
        result = solve_ode_generic(cfg, cond);
        if is_coupled
            p_trace = result.y(18 + n_genes + gi, :);
            plot(ax, result.t, p_trace, 'Color', get_cond_color(cfg, cond), 'LineWidth', 1.5);
        else
            u = compute_u_generic(cfg, result, cond);
            if isfield(u, gene)
                ktx = cfg.ktx_fits.(gene);
                [t_det, ~, p_det] = solve_deterministic_expression(u.t, u.(gene), ktx, ...
                    cfg.expression.beta_m, cfg.expression.beta_p, cfg.expression.ktl, t_eval);
                plot(ax, t_det, p_det, 'Color', get_cond_color(cfg, cond), 'LineWidth', 1.5);
            end
        end
    end
    title(ax, gene, 'FontWeight', 'bold');
    xlim(ax, cfg.solver.t_span);
    xlabel(ax, 'Time (h)'); ylabel(ax, 'Protein (nM)');
    if gi == 1, legend(ax, conditions, 'Location', 'northeast', 'FontSize', 7); end
end
sgtitle(sprintf('%s — Enzyme / Protein Dynamics (coupled ODE)', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig9, fullfile(fig_dir, 'fig9_protein_dynamics.png'), 'Resolution', 300);
saveas(fig9, fullfile(fig_dir, 'fig9_protein_dynamics.svg'));
close(fig9);


%% FIG 10: ENZYME-METABOLITE COUPLING OVERLAY
if is_coupled
    fprintf('  FIG 10: Enzyme-metabolite coupling ...\n');
    % Show key metabolites + their catalyzing enzyme protein in same panel
    % manA->R1(F6P), gmd->R4(gdpman), wbgL->R6(gdpfuc), afcA->R7(fl2_c)
    overlay_genes = [1, 4, 6, 7];  % gene indices
    overlay_met   = [1, 4, 6, 8];  % species indices: F6P, gdpman, gdpfuc, fl2_c
    overlay_labels = {'F6P / ManA', 'GDP-Man / Gmd', 'GDP-Fuc / WbgL', '2-FL / AfcA'};
    fig10 = figure('Position',[100 100 500*2 400*2],'Visible','off');
    for oi = 1:4
        gi_ov = overlay_genes(oi);
        mi_ov = overlay_met(oi);
        for ci = 1:min(n_cond, 4)
            cond = conditions{ci};
            result = solve_ode_generic(cfg, cond);
            ax = subplot(4, min(n_cond,4), (oi-1)*min(n_cond,4) + ci);
            yyaxis(ax, 'left');
            plot(ax, result.t, result.y(mi_ov,:), 'Color', get_cond_color(cfg, cond), 'LineWidth', 2);
            ylabel(ax, 'Metabolite (mM)');
            yyaxis(ax, 'right');
            plot(ax, result.t, result.y(18 + n_genes + gi_ov,:), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
            ylabel(ax, 'Protein (nM)');
            xlim(ax, cfg.solver.t_span);
            if oi == 1, title(ax, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond)); end
            if ci == 1, text(ax, -0.35, 0.5, overlay_labels{oi}, 'Units', 'normalized', 'Rotation', 90, 'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'FontSize', 9); end
            if oi == 4, xlabel(ax, 'Time (h)'); end
        end
    end
    sgtitle(sprintf('%s — Enzyme-Metabolite Coupling', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
    exportgraphics(fig10, fullfile(fig_dir, 'fig10_enzyme_metabolite_coupling.png'), 'Resolution', 300);
    saveas(fig10, fullfile(fig_dir, 'fig10_enzyme_metabolite_coupling.svg'));
    close(fig10);
end

%% FIG 11: S MATRIX HEATMAP
try
fprintf('  FIG 11: Stoichiometric matrix heatmap ...\n');
fig11h = figure('Position',[100 100 450*n_cond 400],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    model = cfg.models.(cond);
    S = model.S;
    ax = subplot(1, n_cond, ci);
    imagesc(ax, S);
    colormap(ax, [0.8 0.2 0.2; 1 1 1; 0.2 0.5 0.8]);
    caxis(ax, [-max(abs(S(:))) max(abs(S(:)))]);
    set(ax, 'YTick', 1:numel(model.species), 'YTickLabel', strrep(model.species, '_', '\_'), 'FontSize', 6);
    set(ax, 'XTick', 1:size(S,2), 'XTickLabel', arrayfun(@(x) sprintf('v%d',x), 1:size(S,2), 'UniformOutput', false), 'FontSize', 6);
    xtickangle(ax, 45);
    title(ax, sprintf('%s (%dx%d)', get_cond_label(cfg, cond), size(S,1), size(S,2)), ...
        'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    if ci == n_cond, colorbar(ax); end
end
sgtitle(sprintf('%s — Stoichiometric Matrices', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig11h, fullfile(fig_dir, 'fig11_S_matrix_heatmap.png'), 'Resolution', 300);
saveas(fig11h, fullfile(fig_dir, 'fig11_S_matrix_heatmap.svg'));
close(fig11h);
catch ME11h
    fprintf('  WARN FIG 11: %s\n', ME11h.message);
end

%% FIG 12: METABOLIC INTERMEDIATES
try
fprintf('  FIG 12: Metabolic intermediates ...\n');
fig12i = figure('Position',[100 100 500*n_cond 500],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    model = cfg.models.(cond);
    result = solve_ode_generic(cfg, cond);
    t = result.t;
    n_sp = numel(model.species);
    skip_idx = [model.substrate_idx, model.target_idx, model.biomass_idx];
    int_idx = setdiff(1:n_sp, skip_idx);
    ax = subplot(1, n_cond, ci);
    hold(ax, 'on');
    cmap = parula(numel(int_idx));
    leg_names = {};
    for ii = 1:numel(int_idx)
        si = int_idx(ii);
        y_sp = result.y(si,:);
        if max(y_sp) > 1e-6
            plot(ax, t, y_sp, 'Color', cmap(ii,:), 'LineWidth', 1.2);
            leg_names{end+1} = strrep(model.species{si}, '_', '\_');
        end
    end
    xlabel(ax, 'Time (h)'); ylabel(ax, 'Concentration (mM)');
    title(ax, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    xlim(ax, cfg.solver.t_span);
    if ~isempty(leg_names), legend(ax, leg_names, 'Location', 'northeast', 'FontSize', 5); end
end
sgtitle(sprintf('%s — Metabolic Intermediates', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig12i, fullfile(fig_dir, 'fig12_intermediates.png'), 'Resolution', 300);
saveas(fig12i, fullfile(fig_dir, 'fig12_intermediates.svg'));
close(fig12i);
catch ME12i
    fprintf('  WARN FIG 12: %s\n', ME12i.message);
end

%% FIG 13: FLUX TIME PROFILES
try
fprintf('  FIG 13: Flux profiles ...\n');
fig13f = figure('Position',[100 100 500*n_cond 500],'Visible','off');
for ci = 1:n_cond
    cond = conditions{ci};
    model = cfg.models.(cond);
    result = solve_ode_generic(cfg, cond);
    t = result.t;
    n_rxn = size(model.S, 2);
    n_t = numel(t);
    flux_mat = zeros(n_rxn, n_t);
    for ti = 1:n_t
        x_t = result.y(:, ti);
        v_all = model.rate_fn(x_t, model.params);
        if isstruct(v_all), v_all = v_all.v; end
        flux_mat(:, ti) = v_all(:);
    end
    ax = subplot(1, n_cond, ci);
    hold(ax, 'on');
    cmap = turbo(n_rxn);
    leg_names = {};
    for ri = 1:n_rxn
        plot(ax, t, flux_mat(ri,:), 'Color', cmap(ri,:), 'LineWidth', 1.2);
        leg_names{end+1} = sprintf('v%d', ri);
    end
    xlabel(ax, 'Time (h)'); ylabel(ax, 'Flux (mM/h)');
    title(ax, get_cond_label(cfg, cond), 'FontWeight', 'bold', 'Color', get_cond_color(cfg, cond));
    xlim(ax, cfg.solver.t_span);
    legend(ax, leg_names, 'Location', 'northeast', 'FontSize', 5);
end
sgtitle(sprintf('%s — Reaction Flux Profiles', cfg.organism), 'FontSize', 12, 'FontWeight', 'bold');
exportgraphics(fig13f, fullfile(fig_dir, 'fig13_flux_profiles.png'), 'Resolution', 300);
saveas(fig13f, fullfile(fig_dir, 'fig13_flux_profiles.svg'));
close(fig13f);
catch ME13f
    fprintf('  WARN FIG 13: %s\n', ME13f.message);
end

n_figs = numel(dir(fullfile(fig_dir, '*.png')));
fprintf('  Generated %d figures\n', n_figs);
end

%% ═══ HELPERS ═══
function c = get_cond_color(cfg, cond)
    if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'colors') && isfield(cfg.visualization.colors, cond)
        c = cfg.visualization.colors.(cond);
    else
        colors = lines(numel(cfg.conditions));
        idx = find(strcmp(cfg.conditions, cond));
        if ~isempty(idx), c = colors(idx,:); else, c = [0 0 0]; end
    end
end

function l = get_cond_label(cfg, cond)
    if isfield(cfg, 'visualization') && isfield(cfg.visualization, 'labels') && isfield(cfg.visualization.labels, cond)
        l = cfg.visualization.labels.(cond);
    else
        l = cond;
    end
end