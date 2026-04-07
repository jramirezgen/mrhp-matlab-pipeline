function results = decipher_routes(gem, target_cpd, opts)
%DECIPHER_ROUTES Discover metabolic routes via BFS/DFS on a GEM network.
%
%   results = DECIPHER_ROUTES(gem, target_cpd, opts)
%
%   Builds a bipartite graph (metabolites <-> reactions) from the GEM
%   stoichiometric matrix, then searches for connected routes between seed
%   metabolites and the target compound. Generates multiple ranked
%   hypotheses with complete route_map, species_map, and S_matrix TSVs.
%
%   Inputs:
%     gem         — struct from load_gem_tsv() with fields S, rxn_ids,
%                   met_ids, rxn_genes, rxn_ec, rxn_equations, etc.
%     target_cpd  — string, target metabolite ID (e.g. 'cpd_MO', 'cpd00020_c0')
%     opts        — struct with optional fields:
%       .mode          — 'degradation' | 'production' (default: 'degradation')
%       .seeds         — cell array of seed metabolite IDs
%                        (default: central metabolism metabolites)
%       .max_hypotheses — max number of routes to return (default: 5)
%       .max_depth      — max BFS depth (default: 15)
%       .output_dir     — directory to write output files
%       .organism       — organism name for labeling (default: 'Unknown')
%       .condition      — condition key (default: 'COND1')
%
%   Output struct:
%     results.hypotheses    — cell array of hypothesis structs
%     results.n_hypotheses  — number found
%     results.ranking       — table with scores
%     results.network_graph — full explored graph struct
%     results.gap_analysis  — per-hypothesis gap info
%
%   Each hypothesis struct contains:
%     .id, .label, .score_prior
%     .rxn_idx, .rxn_ids, .rxn_names, .rxn_equations, .rxn_genes, .rxn_ec
%     .met_idx, .met_ids
%     .S_sub        — sub-network stoichiometric matrix
%     .route_map    — table (8 columns: ode_rxn_idx..note)
%     .species_map  — table (6 columns: ode_idx..note)
%     .is_connected — boolean from gap validation
%     .gaps         — gap struct array
%
%   Example:
%     gem = load_gem_tsv('inputs/shewanella/gem_tsv/');
%     opts.mode = 'degradation';
%     opts.seeds = {'cpd00020_c0','cpd00022_c0'};  % Pyruvate, Succinate
%     opts.output_dir = 'outputs/decipher_MO';
%     results = decipher_routes(gem, 'cpd_MO', opts);
%
%   See also: LOAD_GEM_TSV, VALIDATE_ROUTE_GAPS, BUILD_HYPOTHESIS_CONFIG

    % =====================================================================
    %  PARSE OPTIONS
    % =====================================================================
    if nargin < 3, opts = struct(); end

    mode           = get_opt(opts, 'mode', 'degradation');
    max_hyp        = get_opt(opts, 'max_hypotheses', 5);
    max_depth      = get_opt(opts, 'max_depth', 15);
    output_dir     = get_opt(opts, 'output_dir', '');
    organism_name  = get_opt(opts, 'organism', 'Unknown');
    cond_key       = get_opt(opts, 'condition', 'COND1');

    % Default seeds = central metabolism metabolites
    default_seeds = {'cpd00020_c0', 'cpd00022_c0', 'cpd00027_c0', ...  % Pyr, Succinate, Glucose
                     'cpd00024_c0', 'cpd00032_c0', ...                  % Acetaldehyde, Oxaloacetate
                     'cpd00010_c0', 'cpd00009_c0'};                     % CoA, Pi
    seed_ids = get_opt(opts, 'seeds', default_seeds);

    fprintf('\n=== DECIPHER ROUTES ===\n');
    fprintf('  Target: %s\n', target_cpd);
    fprintf('  Mode:   %s\n', mode);
    fprintf('  Seeds:  %s\n', strjoin(seed_ids, ', '));
    fprintf('  GEM:    %d mets × %d rxns\n', gem.n_mets, gem.n_rxns);

    % =====================================================================
    %  LOCATE TARGET AND SEEDS IN GEM
    % =====================================================================
    target_gem_idx = find_met_idx(gem, target_cpd);
    if isempty(target_gem_idx)
        error('MRHP:decipher:TargetNotFound', ...
              'Target metabolite ''%s'' not found in GEM. Available: %s', ...
              target_cpd, strjoin(gem.met_ids(1:min(10,end)), ', '));
    end

    seed_gem_idx = [];
    valid_seeds = {};
    for i = 1:numel(seed_ids)
        idx = find_met_idx(gem, seed_ids{i});
        if ~isempty(idx)
            seed_gem_idx(end+1) = idx; %#ok<AGROW>
            valid_seeds{end+1}  = seed_ids{i}; %#ok<AGROW>
        end
    end
    fprintf('  Seeds found in GEM: %d / %d\n', numel(seed_gem_idx), numel(seed_ids));

    % =====================================================================
    %  BUILD BIPARTITE ADJACENCY
    % =====================================================================
    % For each reaction, record which metabolites are substrates/products
    fprintf('  Building bipartite graph...\n');

    % met_to_rxn{i} = reactions where metabolite i participates
    met_to_rxn_sub = cell(gem.n_mets, 1);  % as substrate
    met_to_rxn_prod = cell(gem.n_mets, 1); % as product
    for r = 1:gem.n_rxns
        col = gem.S(:, r);
        if issparse(col), col = full(col); end
        subs = find(col < 0);
        prods = find(col > 0);
        for s = subs(:)'
            met_to_rxn_sub{s}(end+1) = r;
        end
        for p = prods(:)'
            met_to_rxn_prod{p}(end+1) = r;
        end
    end

    % =====================================================================
    %  BFS ROUTE SEARCH
    % =====================================================================
    fprintf('  Running BFS route search (max_depth=%d)...\n', max_depth);

    all_paths = {};  % each path = struct with rxn_sequence, met_sequence

    if strcmp(mode, 'degradation')
        % Forward: from target toward central metabolism / seeds
        % Target is consumed → follow products of reactions that consume target
        all_paths = bfs_routes(gem, target_gem_idx, seed_gem_idx, ...
                               met_to_rxn_sub, met_to_rxn_prod, ...
                               max_depth, max_hyp * 3, 'forward');
    else
        % Production: from seeds toward target
        % Seeds are consumed → follow products toward target
        all_paths = bfs_routes(gem, seed_gem_idx, target_gem_idx, ...
                               met_to_rxn_sub, met_to_rxn_prod, ...
                               max_depth, max_hyp * 3, 'forward');
    end

    fprintf('  Found %d candidate routes\n', numel(all_paths));

    if isempty(all_paths)
        warning('MRHP:decipher:NoRoutes', ...
                'No routes found from seeds to target within depth %d.', max_depth);
        results.hypotheses   = {};
        results.n_hypotheses = 0;
        results.ranking      = table();
        results.network_graph = struct();
        results.gap_analysis = {};
        return;
    end

    % =====================================================================
    %  BUILD HYPOTHESES FROM CANDIDATE ROUTES
    % =====================================================================
    fprintf('  Building hypotheses...\n');

    % Score and sort paths
    path_scores = zeros(numel(all_paths), 1);
    for i = 1:numel(all_paths)
        p = all_paths{i};
        % Score = gene coverage × (1 / length) — parsimony + genomic support
        n_rxns_path = numel(p.rxn_sequence);
        n_genes = 0;
        for r = p.rxn_sequence(:)'
            if ~isempty(gem.rxn_genes{r}) && ~strcmp(gem.rxn_genes{r}, '')
                n_genes = n_genes + 1;
            end
        end
        gene_frac = n_genes / max(n_rxns_path, 1);
        parsimony = 1 / (1 + n_rxns_path);
        path_scores(i) = 0.6 * gene_frac + 0.4 * parsimony;
    end

    [~, sort_idx] = sort(path_scores, 'descend');
    n_hyp = min(max_hyp, numel(all_paths));

    hypotheses = cell(1, n_hyp);
    gap_analysis = cell(1, n_hyp);

    for h = 1:n_hyp
        pidx = sort_idx(h);
        p = all_paths{pidx};

        % Collect unique reactions and metabolites
        rxn_set = unique(p.rxn_sequence, 'stable');
        met_set = unique(p.met_sequence, 'stable');

        % Build S sub-matrix
        S_sub = full(gem.S(met_set, rxn_set));

        % Build route_map table (8 columns)
        n_r = numel(rxn_set);
        route_map = table();
        route_map.ode_rxn_idx  = (1:n_r)';
        route_map.ode_rxn_name = cell(n_r, 1);
        route_map.equation     = cell(n_r, 1);
        route_map.gem_rxn_id   = cell(n_r, 1);
        route_map.gene         = cell(n_r, 1);
        route_map.ec           = cell(n_r, 1);
        route_map.evidence     = cell(n_r, 1);
        route_map.note         = cell(n_r, 1);

        for ri = 1:n_r
            gr = rxn_set(ri);
            route_map.ode_rxn_name{ri} = sprintf('v%d_%s', ri, sanitize_id(gem.rxn_ids{gr}));
            route_map.equation{ri}     = gem.rxn_equations{gr};
            route_map.gem_rxn_id{ri}   = gem.rxn_ids{gr};
            route_map.gene{ri}         = gem.rxn_genes{gr};
            route_map.ec{ri}           = gem.rxn_ec{gr};
            if ~isempty(gem.rxn_genes{gr}) && ~strcmp(gem.rxn_genes{gr}, '')
                route_map.evidence{ri} = 'PRESENT_IN_GEM';
            else
                route_map.evidence{ri} = 'INFERRED';
            end
            route_map.note{ri} = 'auto-deciphered';
        end

        % Build species_map table (6 columns)
        n_m = numel(met_set);
        species_map = table();
        species_map.ode_idx      = (1:n_m)';
        species_map.ode_name     = cell(n_m, 1);
        species_map.gem_cpd_id   = cell(n_m, 1);
        species_map.formula      = cell(n_m, 1);
        species_map.compartment  = cell(n_m, 1);
        species_map.note         = cell(n_m, 1);

        for mi = 1:n_m
            gm = met_set(mi);
            species_map.ode_name{mi}    = sanitize_id(gem.met_ids{gm});
            species_map.gem_cpd_id{mi}  = gem.met_ids{gm};
            if gm <= numel(gem.met_formulas)
                species_map.formula{mi} = gem.met_formulas{gm};
            else
                species_map.formula{mi} = '';
            end
            if gm <= numel(gem.met_compartments)
                species_map.compartment{mi} = gem.met_compartments{gm};
            else
                species_map.compartment{mi} = 'c0';
            end
            species_map.note{mi} = 'auto-deciphered';
        end

        % Validate gaps
        % Map target and seed to sub-network indices
        target_sub = find(met_set == target_gem_idx, 1);
        seed_sub   = [];
        for si = 1:numel(seed_gem_idx)
            idx_s = find(met_set == seed_gem_idx(si), 1);
            if ~isempty(idx_s)
                seed_sub(end+1) = idx_s; %#ok<AGROW>
            end
        end
        if isempty(seed_sub), seed_sub = 1; end
        if isempty(target_sub), target_sub = n_m; end

        [is_conn, gaps, gap_report] = validate_route_gaps(S_sub, ...
            species_map.gem_cpd_id, route_map.gem_rxn_id, target_sub, seed_sub);

        % Collect diagnostic genes
        diag_genes = {};
        for ri = 1:n_r
            g = gem.rxn_genes{rxn_set(ri)};
            if ~isempty(g) && ~strcmp(g, '')
                % Split gene_rule by 'and'/'or' and clean
                parts = regexp(g, '\s+(and|or)\s+', 'split');
                for pp = 1:numel(parts)
                    gname = strtrim(parts{pp});
                    gname = regexprep(gname, '[()]', '');
                    if ~isempty(gname) && ~any(strcmp(diag_genes, gname))
                        diag_genes{end+1} = gname; %#ok<AGROW>
                    end
                end
            end
        end

        % Build hypothesis struct
        hyp = struct();
        hyp.id             = sprintf('HYP_%s_%d', cond_key, h);
        hyp.label          = sprintf('Route %d (%d rxns, %d genes)', h, n_r, numel(diag_genes));
        hyp.score_prior    = path_scores(pidx);
        hyp.rxn_idx        = rxn_set;
        hyp.rxn_ids        = gem.rxn_ids(rxn_set);
        hyp.rxn_names      = route_map.ode_rxn_name;
        hyp.rxn_equations  = gem.rxn_equations(rxn_set);
        hyp.rxn_genes      = gem.rxn_genes(rxn_set);
        hyp.rxn_ec         = gem.rxn_ec(rxn_set);
        hyp.met_idx        = met_set;
        hyp.met_ids        = gem.met_ids(met_set);
        hyp.S_sub          = S_sub;
        hyp.route_map      = route_map;
        hyp.species_map    = species_map;
        hyp.is_connected   = is_conn;
        hyp.gaps           = gaps;
        hyp.gap_report     = gap_report;
        hyp.diagnostic_up  = diag_genes;
        hyp.n_reactions    = n_r;
        hyp.n_metabolites  = n_m;
        hyp.n_genes        = numel(diag_genes);

        hypotheses{h} = hyp;
        gap_analysis{h} = struct('hypothesis', h, 'is_connected', is_conn, ...
                                 'n_gaps', numel(gaps), 'gaps', gaps);
    end

    % =====================================================================
    %  BUILD RANKING SUMMARY
    % =====================================================================
    ranking = table();
    ranking.hypothesis_id  = cell(n_hyp, 1);
    ranking.label          = cell(n_hyp, 1);
    ranking.n_reactions    = zeros(n_hyp, 1);
    ranking.n_genes        = zeros(n_hyp, 1);
    ranking.n_metabolites  = zeros(n_hyp, 1);
    ranking.gene_coverage  = zeros(n_hyp, 1);
    ranking.parsimony      = zeros(n_hyp, 1);
    ranking.score_prior    = zeros(n_hyp, 1);
    ranking.is_connected   = false(n_hyp, 1);
    ranking.n_gaps         = zeros(n_hyp, 1);

    for h = 1:n_hyp
        hy = hypotheses{h};
        ranking.hypothesis_id{h} = hy.id;
        ranking.label{h}         = hy.label;
        ranking.n_reactions(h)   = hy.n_reactions;
        ranking.n_genes(h)       = hy.n_genes;
        ranking.n_metabolites(h) = hy.n_metabolites;
        ranking.gene_coverage(h) = hy.n_genes / max(hy.n_reactions, 1);
        ranking.parsimony(h)     = 1 / (1 + hy.n_reactions);
        ranking.score_prior(h)   = hy.score_prior;
        ranking.is_connected(h)  = hy.is_connected;
        ranking.n_gaps(h)        = numel(hy.gaps);
    end

    % =====================================================================
    %  EXPORT FILES
    % =====================================================================
    if ~isempty(output_dir)
        fprintf('  Writing output to %s\n', output_dir);

        % Create directories
        if ~isfolder(output_dir), mkdir(output_dir); end
        hyp_root = fullfile(output_dir, 'hypotheses');
        if ~isfolder(hyp_root), mkdir(hyp_root); end
        snap_dir = fullfile(output_dir, 'gem_snapshot');
        if ~isfolder(snap_dir), mkdir(snap_dir); end

        % Write per-hypothesis files
        for h = 1:n_hyp
            hy = hypotheses{h};
            hdir = fullfile(hyp_root, sprintf('hypothesis_%d', h));
            if ~isfolder(hdir), mkdir(hdir); end

            % route_map.tsv
            writetable(hy.route_map, fullfile(hdir, 'route_map.tsv'), ...
                       'FileType','text', 'Delimiter','\t');

            % species_map.tsv
            writetable(hy.species_map, fullfile(hdir, 'species_map.tsv'), ...
                       'FileType','text', 'Delimiter','\t');

            % S_matrix.tsv
            write_s_matrix(hy.S_sub, hy.species_map.gem_cpd_id, ...
                           hy.route_map.ode_rxn_name, ...
                           fullfile(hdir, 'S_matrix.tsv'));

            % rate_fn_template.m
            write_rate_fn_template(hy, fullfile(hdir, 'rate_fn_template.m'));

            % config_hypothesis.m
            write_hypothesis_config(hy, fullfile(hdir, 'config_hypothesis.m'), ...
                                    organism_name, cond_key, opts);
        end

        % ranking_summary.tsv
        writetable(ranking, fullfile(output_dir, 'ranking_summary.tsv'), ...
                   'FileType','text', 'Delimiter','\t');

        % gap_analysis.tsv
        write_gap_analysis(gap_analysis, fullfile(output_dir, 'gap_analysis.tsv'));

        % decipher_report.md
        write_decipher_report(output_dir, target_cpd, mode, hypotheses, ...
                              ranking, organism_name, cond_key);

        fprintf('  Output written: %d hypotheses\n', n_hyp);
    end

    % =====================================================================
    %  ASSEMBLE OUTPUT
    % =====================================================================
    results.hypotheses    = hypotheses;
    results.n_hypotheses  = n_hyp;
    results.ranking       = ranking;
    results.network_graph = struct('target', target_cpd, 'seeds', {valid_seeds}, ...
                                  'mode', mode, 'n_paths_explored', numel(all_paths));
    results.gap_analysis  = gap_analysis;

    fprintf('=== DECIPHER COMPLETE: %d hypotheses generated ===\n\n', n_hyp);
end

% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function paths = bfs_routes(gem, start_mets, end_mets, met_to_rxn_sub, met_to_rxn_prod, max_depth, max_paths, ~)
%BFS_ROUTES BFS from start metabolites, collect paths reaching end metabolites.
    paths = {};

    % Queue entries: struct with current_met, rxn_sequence, met_sequence, visited_mets
    queue = {};
    if isscalar(start_mets)
        start_mets_vec = start_mets;
    else
        start_mets_vec = start_mets(:)';
    end

    for sm = start_mets_vec
        entry.current_met   = sm;
        entry.rxn_sequence  = [];
        entry.met_sequence  = sm;
        entry.visited_mets  = sm;
        entry.visited_rxns  = [];
        queue{end+1} = entry; %#ok<AGROW>
    end

    if isscalar(end_mets)
        end_set = end_mets;
    else
        end_set = end_mets(:)';
    end

    head = 1;
    while head <= numel(queue) && numel(paths) < max_paths
        node = queue{head};
        head = head + 1;

        if numel(node.rxn_sequence) >= max_depth
            continue;
        end

        cm = node.current_met;

        % Find reactions that consume this metabolite
        rxns_consuming = met_to_rxn_sub{cm};

        for r = rxns_consuming(:)'
            if any(node.visited_rxns == r), continue; end

            % Find products of this reaction
            col = gem.S(:, r);
            if issparse(col), col = full(col); end
            prods = find(col > 0);

            for pm = prods(:)'
                if any(node.visited_mets == pm), continue; end  % avoid cycles

                new_entry.current_met  = pm;
                new_entry.rxn_sequence = [node.rxn_sequence, r];
                new_entry.met_sequence = [node.met_sequence, pm];
                new_entry.visited_mets = [node.visited_mets, pm];
                new_entry.visited_rxns = [node.visited_rxns, r];

                % Check if we reached an end metabolite
                if any(end_set == pm)
                    path.rxn_sequence = new_entry.rxn_sequence;
                    path.met_sequence = new_entry.met_sequence;
                    paths{end+1} = path; %#ok<AGROW>
                    if numel(paths) >= max_paths, return; end
                end

                % Continue BFS
                if numel(new_entry.rxn_sequence) < max_depth
                    queue{end+1} = new_entry; %#ok<AGROW>
                end
            end
        end
    end
end

function idx = find_met_idx(gem, met_id)
%FIND_MET_IDX Find metabolite index by exact or partial match.
    idx = find(strcmp(gem.met_ids, met_id), 1);
    if isempty(idx)
        % Try partial match
        matches = find(contains(gem.met_ids, met_id));
        if numel(matches) == 1
            idx = matches;
        elseif numel(matches) > 1
            idx = matches(1);
            fprintf('  Warning: ''%s'' matches %d metabolites, using first: %s\n', ...
                    met_id, numel(matches), gem.met_ids{idx});
        end
    end
end

function s = sanitize_id(id)
%SANITIZE_ID Make an ID safe for MATLAB variable names.
    s = regexprep(id, '[^a-zA-Z0-9_]', '_');
    if ~isempty(s) && ~isstrprop(s(1), 'alpha')
        s = ['x_' s];
    end
end

function v = get_opt(opts, field, default)
%GET_OPT Get optional field or return default.
    if isfield(opts, field)
        v = opts.(field);
    else
        v = default;
    end
end

function write_s_matrix(S, met_ids, rxn_names, filepath)
%WRITE_S_MATRIX Write stoichiometric matrix as TSV.
    fid = fopen(filepath, 'w');
    % Header
    fprintf(fid, 'species\\reaction');
    for j = 1:numel(rxn_names)
        fprintf(fid, '\t%s', rxn_names{j});
    end
    fprintf(fid, '\n');
    % Rows
    for i = 1:size(S, 1)
        fprintf(fid, '%s', met_ids{i});
        for j = 1:size(S, 2)
            fprintf(fid, '\t%g', S(i, j));
        end
        fprintf(fid, '\n');
    end
    fclose(fid);
end

function write_rate_fn_template(hyp, filepath)
%WRITE_RATE_FN_TEMPLATE Write a default Michaelis-Menten rate_fn .m file.
    fid = fopen(filepath, 'w');
    fprintf(fid, 'function v = rate_fn_%s(x, params)\n', lower(hyp.id));
    fprintf(fid, '%%RATE_FN_%s Auto-generated rate function template.\n', upper(hyp.id));
    fprintf(fid, '%%   v = rate_fn_%s(x, params)\n', lower(hyp.id));
    fprintf(fid, '%%\n');
    fprintf(fid, '%%   Default: Michaelis-Menten kinetics for each reaction.\n');
    fprintf(fid, '%%   EDIT THIS FILE to customize kinetics for your system.\n');
    fprintf(fid, '%%\n');
    fprintf(fid, '%%   Species (x indices):\n');
    for mi = 1:numel(hyp.met_ids)
        fprintf(fid, '%%     x(%d) = %s\n', mi, hyp.met_ids{mi});
    end
    fprintf(fid, '%%\n');
    fprintf(fid, '%%   Reactions (v indices):\n');
    for ri = 1:numel(hyp.rxn_ids)
        fprintf(fid, '%%     v(%d) = %s  [%s]\n', ri, hyp.rxn_ids{ri}, hyp.rxn_equations{ri});
    end
    fprintf(fid, '\n');
    fprintf(fid, '    n_rxn = %d;\n', hyp.n_reactions);
    fprintf(fid, '    v = zeros(n_rxn, 1);\n\n');
    fprintf(fid, '    %% Default Michaelis-Menten: v = Vmax * S / (Km + S)\n');
    fprintf(fid, '    Vmax_default = 1.0;\n');
    fprintf(fid, '    Km_default   = 0.5;\n\n');
    for ri = 1:hyp.n_reactions
        col = hyp.S_sub(:, ri);
        subs_idx = find(col < 0);
        if ~isempty(subs_idx)
            si = subs_idx(1);
            fprintf(fid, '    v(%d) = Vmax_default * x(%d) / (Km_default + x(%d));  %% %s\n', ...
                    ri, si, si, hyp.rxn_ids{ri});
        else
            fprintf(fid, '    v(%d) = Vmax_default;  %% %s (no substrate)\n', ri, hyp.rxn_ids{ri});
        end
    end
    fprintf(fid, '\nend\n');
    fclose(fid);
end

function write_hypothesis_config(hyp, filepath, organism, cond_key, opts)
%WRITE_HYPOTHESIS_CONFIG Write a validate-ready config .m file.
    fid = fopen(filepath, 'w');
    fn_name = sprintf('config_hypothesis_%d', str2double(regexp(hyp.id, '\d+$', 'match', 'once')));
    if isnan(str2double(regexp(hyp.id, '\d+$', 'match', 'once')))
        fn_name = 'config_hypothesis';
    end

    fprintf(fid, 'function cfg = %s()\n', fn_name);
    fprintf(fid, '%%%s Configuration auto-generated by MRHP decipher mode.\n', upper(fn_name));
    fprintf(fid, '%%   Generated from hypothesis: %s\n', hyp.id);
    fprintf(fid, '%%   Target: %s | Mode: %s\n', ...
            get_opt(opts, 'target', 'unknown'), get_opt(opts, 'mode', 'degradation'));
    fprintf(fid, '\n');
    fprintf(fid, 'cfg.organism = ''%s (Deciphered)'';\n', organism);
    fprintf(fid, 'cfg.conditions = {''%s''};\n', cond_key);
    fprintf(fid, '\n');

    % Genes
    if ~isempty(hyp.diagnostic_up)
        fprintf(fid, 'cfg.genes = {');
        for gi = 1:numel(hyp.diagnostic_up)
            if gi > 1, fprintf(fid, ', '); end
            fprintf(fid, '''%s''', hyp.diagnostic_up{gi});
        end
        fprintf(fid, '};\n\n');
    else
        fprintf(fid, 'cfg.genes = {''Gene1''};\n\n');
    end

    % Solver defaults
    fprintf(fid, '%%%% SOLVER\n');
    fprintf(fid, 'cfg.solver.t_span = [0, 72];\n');
    fprintf(fid, 'cfg.solver.t_eval = (0:0.5:72)'';\n');
    fprintf(fid, 'cfg.solver.rtol = 1e-8;\n');
    fprintf(fid, 'cfg.solver.atol = 1e-10;\n');
    fprintf(fid, 'cfg.solver.max_step = 0.5;\n');
    fprintf(fid, 'cfg.solver.substrate_conversion = 1.0;\n\n');

    % Sweep defaults
    fprintf(fid, '%%%% SWEEP\n');
    fprintf(fid, 'cfg.sweep.substrate_values = [1.0];\n');
    fprintf(fid, 'cfg.sweep.target_values = [1.0];\n');
    fprintf(fid, 'cfg.ref_substrate = 1.0;\n');
    fprintf(fid, 'cfg.ref_target = 1.0;\n');
    fprintf(fid, 'cfg.t_sample = 24.0;\n\n');

    % Expression defaults
    fprintf(fid, '%%%% EXPRESSION\n');
    fprintf(fid, 'cfg.expression.beta_m = 8.32;\n');
    fprintf(fid, 'cfg.expression.beta_p = 0.924;\n');
    fprintf(fid, 'cfg.expression.ktl = 6.0;\n');
    fprintf(fid, 'cfg.expression.stochastic_tau = 0.01;\n');
    fprintf(fid, 'cfg.expression.n_cells = 20;\n');
    fprintf(fid, 'cfg.expression.nm_to_mol = 6022;\n\n');

    % Ktx defaults
    fprintf(fid, '%%%% TRANSCRIPTION RATES (edit with experimental values)\n');
    for gi = 1:numel(hyp.diagnostic_up)
        fprintf(fid, 'cfg.ktx_fits.%s = 0.1;\n', hyp.diagnostic_up{gi});
    end
    fprintf(fid, '\n');

    % Model: species, S, x0
    fprintf(fid, '%%%% MODEL — %s\n', cond_key);
    fprintf(fid, 'm.species = {');
    for mi = 1:hyp.n_metabolites
        if mi > 1, fprintf(fid, ', '); end
        fprintf(fid, '''%s''', hyp.species_map.ode_name{mi});
    end
    fprintf(fid, '};\n\n');

    fprintf(fid, 'm.S = [...\n');
    for i = 1:size(hyp.S_sub, 1)
        fprintf(fid, '  ');
        for j = 1:size(hyp.S_sub, 2)
            fprintf(fid, '%3g', hyp.S_sub(i, j));
            if j < size(hyp.S_sub, 2), fprintf(fid, ','); end
        end
        if i < size(hyp.S_sub, 1)
            fprintf(fid, ';...\n');
        else
            fprintf(fid, '];\n\n');
        end
    end

    % Initial conditions (1.0 for first species, 0.1 for rest)
    fprintf(fid, 'm.base_x0 = [');
    for mi = 1:hyp.n_metabolites
        if mi == 1
            fprintf(fid, '1.0');
        else
            fprintf(fid, '; 0.1');
        end
    end
    fprintf(fid, '];\n\n');

    % Indices
    fprintf(fid, 'm.substrate_idx = 1;\n');
    fprintf(fid, 'm.target_idx = %d;  %% EDIT: set to actual target species index\n', ...
            min(hyp.n_metabolites, 2));
    fprintf(fid, 'm.biomass_idx = %d;  %% EDIT: set to actual biomass species index\n', ...
            hyp.n_metabolites);
    fprintf(fid, '\n');

    % Rate function reference
    fprintf(fid, 'm.params = struct();\n');
    fprintf(fid, 'm.rate_fn = @rate_fn_%s;  %% See rate_fn_template.m in same directory\n\n', lower(hyp.id));
    fprintf(fid, 'cfg.models.%s = m;\n\n', cond_key);

    % Phi definition
    fprintf(fid, '%%%% PHI DEFINITION\n');
    fprintf(fid, 'cfg.phi_definition.default = struct(''type'',''depletion'',''target_idx'',m.target_idx);\n\n');

    % Phenotype reference (placeholder)
    fprintf(fid, '%%%% PHENOTYPE REFERENCE (EDIT with experimental data)\n');
    fprintf(fid, 'cfg.phenotype_ref.%s = struct(''type'',''exponential_depletion'',''mu'',0.1);\n\n', cond_key);

    % Regulatory (placeholder)
    fprintf(fid, '%%%% REGULATORY SIGNALS (EDIT per gene)\n');
    for gi = 1:numel(hyp.diagnostic_up)
        fprintf(fid, 'cfg.regulatory.%s.default = struct(''form'',''act'',''signals'',{{1}},''K'',[0.5],''n'',[2],''basal'',0.01);\n', ...
                hyp.diagnostic_up{gi});
    end
    fprintf(fid, '\n');

    % Hypotheses
    fprintf(fid, '%%%% HYPOTHESES\n');
    fprintf(fid, 'cfg.hypotheses.%s = {struct(''id'',''%s'',''label'',''%s'',...\n', ...
            cond_key, hyp.id, hyp.label);
    fprintf(fid, '    ''diagnostic_up'',{{');
    for gi = 1:numel(hyp.diagnostic_up)
        if gi > 1, fprintf(fid, ','); end
        fprintf(fid, '''%s''', hyp.diagnostic_up{gi});
    end
    fprintf(fid, '}},''score_total'',%.4f)};\n\n', hyp.score_prior);

    % Route maps
    fprintf(fid, '%%%% ROUTE MAPS\n');
    fprintf(fid, 'hyp_dir = fileparts(mfilename(''fullpath''));\n');
    fprintf(fid, 'cfg.route_maps.%s = struct(...\n', cond_key);
    fprintf(fid, '    ''route_map'', fullfile(hyp_dir, ''route_map.tsv''),...\n');
    fprintf(fid, '    ''species_map'', fullfile(hyp_dir, ''species_map.tsv''),...\n');
    fprintf(fid, '    ''S_matrix'', fullfile(hyp_dir, ''S_matrix.tsv''));\n\n');

    % Visualization
    fprintf(fid, '%%%% VISUALIZATION\n');
    fprintf(fid, 'cfg.visualization.colors.%s = [0.2, 0.6, 0.8];\n', cond_key);
    fprintf(fid, 'cfg.visualization.labels.%s = ''%s'';\n\n', cond_key, cond_key);

    fprintf(fid, 'end\n');
    fclose(fid);
end

function write_gap_analysis(gap_analysis, filepath)
%WRITE_GAP_ANALYSIS Write gap analysis TSV.
    fid = fopen(filepath, 'w');
    fprintf(fid, 'hypothesis\tis_connected\tn_gaps\tgap_metabolite\tgap_type\n');
    for h = 1:numel(gap_analysis)
        ga = gap_analysis{h};
        if ga.n_gaps == 0
            fprintf(fid, '%d\t%d\t%d\t-\t-\n', h, ga.is_connected, 0);
        else
            for g = 1:numel(ga.gaps)
                fprintf(fid, '%d\t%d\t%d\t%s\t%s\n', h, ga.is_connected, ga.n_gaps, ...
                        ga.gaps(g).metabolite, ga.gaps(g).type);
            end
        end
    end
    fclose(fid);
end

function write_decipher_report(output_dir, target_cpd, mode, hypotheses, ranking, organism, cond_key)
%WRITE_DECIPHER_REPORT Write markdown report.
    fid = fopen(fullfile(output_dir, 'decipher_report.md'), 'w');
    fprintf(fid, '# MRHP Decipher Report\n\n');
    fprintf(fid, '**Generated:** %s\n\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '| Field | Value |\n|-------|-------|\n');
    fprintf(fid, '| Organism | %s |\n', organism);
    fprintf(fid, '| Condition | %s |\n', cond_key);
    fprintf(fid, '| Target | `%s` |\n', target_cpd);
    fprintf(fid, '| Mode | %s |\n', mode);
    fprintf(fid, '| Hypotheses | %d |\n\n', numel(hypotheses));

    fprintf(fid, '## Ranking Summary\n\n');
    fprintf(fid, '| Rank | ID | Reactions | Genes | Connected | Score |\n');
    fprintf(fid, '|------|----|-----------|-------|-----------|-------|\n');
    for h = 1:height(ranking)
        fprintf(fid, '| %d | %s | %d | %d | %s | %.4f |\n', ...
                h, ranking.hypothesis_id{h}, ranking.n_reactions(h), ...
                ranking.n_genes(h), mat2str(ranking.is_connected(h)), ...
                ranking.score_prior(h));
    end
    fprintf(fid, '\n');

    for h = 1:numel(hypotheses)
        hy = hypotheses{h};
        fprintf(fid, '## Hypothesis %d: %s\n\n', h, hy.label);
        fprintf(fid, '**Diagnostic genes:** %s\n\n', strjoin(hy.diagnostic_up, ', '));
        fprintf(fid, '**Reactions:**\n\n');
        for ri = 1:hy.n_reactions
            fprintf(fid, '%d. `%s` — %s\n', ri, hy.rxn_ids{ri}, hy.rxn_equations{ri});
        end
        fprintf(fid, '\n**Gap status:** %s (%d gaps)\n\n', ...
                ternary(hy.is_connected, 'CONNECTED', 'GAPS DETECTED'), numel(hy.gaps));
        if ~isempty(hy.gaps)
            for g = 1:numel(hy.gaps)
                fprintf(fid, '- GAP: `%s` (%s)\n', hy.gaps(g).metabolite, hy.gaps(g).type);
            end
            fprintf(fid, '\n');
        end
        fprintf(fid, '---\n\n');
    end

    fprintf(fid, '## Next Steps\n\n');
    fprintf(fid, '1. Review and edit `config_hypothesis.m` in each hypothesis folder\n');
    fprintf(fid, '2. Customize `rate_fn_template.m` with system-specific kinetics\n');
    fprintf(fid, '3. Run validation: `mrhp(''validate'', ''config'', ''path/to/config_hypothesis.m'', ...)`\n');
    fclose(fid);
end

function r = ternary(cond, a, b)
    if cond, r = a; else, r = b; end
end
