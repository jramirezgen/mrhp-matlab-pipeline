function mrhp(mode, varargin)
%MRHP Multiscale Route Hypothesis Platform — Unified Entry Point.
%
%   MRHP(mode, Name, Value, ...) runs the MRHP pipeline in one of two modes:
%
%   MODES:
%     'decipher'  — Discover metabolic routes from a GEM via BFS/DFS graph
%                   search. Generates 3-5 ranked hypotheses with complete
%                   route_map, species_map, S_matrix TSVs, rate function
%                   templates, and validate-ready config files.
%
%     'validate'  — Validate hypotheses through the full 13-phase pipeline:
%                   ODE simulation, phi extraction, bridge fitting, regulatory
%                   signals, deterministic + stochastic expression, hypothesis
%                   scoring, and 13 publication-quality figures.
%
%     'help'      — Display this help message with usage examples.
%
%   ─── DECIPHER MODE ─────────────────────────────────────────────────────
%
%   mrhp('decipher', Name, Value, ...)
%
%   Required arguments:
%     'gem_dir'      — Path to directory with GEM TSV exports:
%                      reactions.tsv, metabolites.tsv, S_matrix.tsv
%     'target'       — Target metabolite ID (e.g. 'cpd_MO', 'cpd00020_c0')
%     'output_dir'   — Directory for output files
%
%   Optional arguments:
%     'mode'            — 'degradation' (default) | 'production'
%                         Degradation: target → central metabolism
%                         Production:  seeds  → target
%     'seeds'           — Cell array of seed metabolite IDs
%                         Default: central metabolism (Pyr, Succinate, etc.)
%     'max_hypotheses'  — Maximum hypotheses to generate (default: 5)
%     'max_depth'       — Maximum BFS search depth (default: 15)
%     'organism'        — Organism name for labeling (default: 'Unknown')
%     'condition'       — Condition key (default: 'COND1')
%
%   Output structure (in output_dir):
%     hypotheses/hypothesis_N/   — Per-hypothesis files:
%       route_map.tsv              8 cols: ode_rxn_idx, name, equation,
%                                         gem_rxn_id, gene, ec, evidence, note
%       species_map.tsv            6 cols: ode_idx, name, gem_cpd_id,
%                                         formula, compartment, note
%       S_matrix.tsv               Stoichiometric matrix
%       config_hypothesis.m        Ready for validate mode
%       rate_fn_template.m         Michaelis-Menten template
%     ranking_summary.tsv        — Comparative ranking of all hypotheses
%     gap_analysis.tsv           — Gap detection per hypothesis
%     decipher_report.md         — Full markdown report
%
%   ─── VALIDATE MODE ─────────────────────────────────────────────────────
%
%   mrhp('validate', Name, Value, ...)
%
%   Required arguments:
%     'config'       — Configuration source. One of:
%                      - Function handle:   @config_shewanella_lc6
%                      - Path to .m file:   'path/to/config_hypothesis.m'
%                      - Path to hyp dir:   'outputs/decipher/hypotheses/hypothesis_1'
%     'output_dir'   — Directory for pipeline output
%
%   Optional arguments:
%     'experimental_data' — Path to experimental data directory
%                           (for RT-qPCR, kinetics, etc.)
%     'pipeline_mode'     — 'full' (default) | 'quick' | 'figures_only'
%     'organism'          — Override organism name
%     'condition'         — Override condition key
%
%   Output structure (in output_dir):
%     ode_timeseries/    — ODE concentration traces (TSV)
%     phi_signals/       — Metabolic activity signals (TSV)
%     bridge_predictions/— Bridge fit predictions + R² (TSV)
%     regulatory_signals/— Per-gene regulatory signals (TSV)
%     expression/        — Deterministic + stochastic expression (TSV)
%     hypothesis_ranking/— 5-component scores (TSV)
%     route_maps/        — GEM traceability (TSV + equations.md)
%     figures/           — 13 publication figures (PNG 300dpi + SVG)
%     metrics/           — Bridge metrics per scenario
%     summary_metrics/   — Full grid summary
%     metadata/          — run_info.json, certification.json
%     reports/           — pipeline_execution_report.md
%
%   ─── EXAMPLES ──────────────────────────────────────────────────────────
%
%   % Discover routes for Methyl Orange degradation
%   mrhp('decipher', ...
%        'gem_dir',    'inputs/shewanella/gem_tsv/', ...
%        'target',     'cpd_MO', ...
%        'mode',       'degradation', ...
%        'seeds',      {'cpd_YE_eff'}, ...
%        'max_hypotheses', 5, ...
%        'output_dir', 'outputs/decipher_MO')
%
%   % Validate from decipher output
%   mrhp('validate', ...
%        'config', 'outputs/decipher_MO/hypotheses/hypothesis_1', ...
%        'experimental_data', 'inputs/shewanella/', ...
%        'output_dir', 'outputs/validate_MO')
%
%   % Validate from existing config
%   mrhp('validate', ...
%        'config', @config_shewanella_lc6, ...
%        'output_dir', 'outputs/validate_shewanella')
%
%   % Quick validation (skip parameter grid)
%   mrhp('validate', ...
%        'config', @config_ecoli_fucose, ...
%        'output_dir', 'outputs/validate_ecoli', ...
%        'pipeline_mode', 'quick')
%
%   See also: DECIPHER_ROUTES, LOAD_GEM_TSV, VALIDATE_ROUTE_GAPS,
%             BUILD_HYPOTHESIS_CONFIG, RUN_PIPELINE_GENERIC

    % =====================================================================
    %  SETUP
    % =====================================================================
    if nargin == 0, mode = 'help'; end

    % Add engine and configs to path
    base_dir = fileparts(mfilename('fullpath'));
    addpath(fullfile(base_dir, 'engine'));
    addpath(fullfile(base_dir, 'configs'));
    addpath(fullfile(base_dir, 'configs', 'route_maps'));

    fprintf('\n');
    fprintf('╔══════════════════════════════════════════════════════════╗\n');
    fprintf('║  MRHP — Multiscale Route Hypothesis Platform  v5.0.0   ║\n');
    fprintf('╚══════════════════════════════════════════════════════════╝\n\n');

    % =====================================================================
    %  DISPATCH
    % =====================================================================
    switch lower(mode)
        case 'help'
            help(mfilename);
            return;

        case 'decipher'
            run_decipher(varargin);

        case 'validate'
            run_validate(varargin);

        otherwise
            error('MRHP:InvalidMode', ...
                  'Unknown mode ''%s''. Use ''decipher'', ''validate'', or ''help''.', mode);
    end
end

% =========================================================================
%  DECIPHER MODE
% =========================================================================

function run_decipher(args)
    p = inputParser;
    p.FunctionName = 'mrhp(''decipher'', ...)';

    addParameter(p, 'gem_dir',         '',    @ischar);
    addParameter(p, 'target',          '',    @ischar);
    addParameter(p, 'output_dir',      '',    @ischar);
    addParameter(p, 'mode',            'degradation', @ischar);
    addParameter(p, 'seeds',           {},    @iscell);
    addParameter(p, 'max_hypotheses',  5,     @isnumeric);
    addParameter(p, 'max_depth',       15,    @isnumeric);
    addParameter(p, 'organism',        'Unknown', @ischar);
    addParameter(p, 'condition',       'COND1',   @ischar);

    parse(p, args{:});
    r = p.Results;

    % Validate required
    if isempty(r.gem_dir)
        error('MRHP:decipher:NoGemDir', ...
              'Missing required argument ''gem_dir''. Use mrhp(''help'') for usage.');
    end
    if isempty(r.target)
        error('MRHP:decipher:NoTarget', ...
              'Missing required argument ''target''. Use mrhp(''help'') for usage.');
    end
    if isempty(r.output_dir)
        error('MRHP:decipher:NoOutputDir', ...
              'Missing required argument ''output_dir''. Use mrhp(''help'') for usage.');
    end

    fprintf('MODE: DECIPHER\n');
    fprintf('  GEM dir:    %s\n', r.gem_dir);
    fprintf('  Target:     %s\n', r.target);
    fprintf('  Search:     %s (depth=%d, max=%d hypotheses)\n', r.mode, r.max_depth, r.max_hypotheses);
    fprintf('  Output:     %s\n\n', r.output_dir);

    % Load GEM
    gem = load_gem_tsv(r.gem_dir);

    % Build opts for decipher_routes
    opts = struct();
    opts.mode           = r.mode;
    opts.max_hypotheses = r.max_hypotheses;
    opts.max_depth      = r.max_depth;
    opts.output_dir     = r.output_dir;
    opts.organism       = r.organism;
    opts.condition      = r.condition;
    opts.target         = r.target;
    if ~isempty(r.seeds)
        opts.seeds = r.seeds;
    end

    % Run decipher
    results = decipher_routes(gem, r.target, opts);

    % Summary
    fprintf('\n╔══════════════════════════════════════════════════════════╗\n');
    fprintf('║  DECIPHER COMPLETE                                     ║\n');
    fprintf('╠══════════════════════════════════════════════════════════╣\n');
    fprintf('║  Hypotheses generated: %d                               ║\n', results.n_hypotheses);
    for h = 1:results.n_hypotheses
        hy = results.hypotheses{h};
        fprintf('║  %d. %s (score=%.3f, connected=%d)  \n', ...
                h, hy.id, hy.score_prior, hy.is_connected);
    end
    fprintf('║                                                         ║\n');
    fprintf('║  Output: %s\n', results.hypotheses{1}.id);
    fprintf('║  Next: mrhp(''validate'', ''config'', ''<hyp_dir>'', ...)     ║\n');
    fprintf('╚══════════════════════════════════════════════════════════╝\n\n');
end

% =========================================================================
%  VALIDATE MODE
% =========================================================================

function run_validate(args)
    p = inputParser;
    p.FunctionName = 'mrhp(''validate'', ...)';

    addParameter(p, 'config',            [], @(x) true);  % fn handle, char, or string
    addParameter(p, 'output_dir',        '', @ischar);
    addParameter(p, 'experimental_data', '', @ischar);
    addParameter(p, 'pipeline_mode',     'full', @ischar);
    addParameter(p, 'organism',          '', @ischar);
    addParameter(p, 'condition',         '', @ischar);

    parse(p, args{:});
    r = p.Results;

    if isempty(r.config)
        error('MRHP:validate:NoConfig', ...
              'Missing required argument ''config''. Use mrhp(''help'') for usage.');
    end
    if isempty(r.output_dir)
        error('MRHP:validate:NoOutputDir', ...
              'Missing required argument ''output_dir''. Use mrhp(''help'') for usage.');
    end

    fprintf('MODE: VALIDATE\n');

    % --- Resolve config ---
    if isa(r.config, 'function_handle')
        % Direct function handle: @config_shewanella_lc6
        fprintf('  Config: function handle %s\n', func2str(r.config));
        cfg = r.config();

    elseif ischar(r.config) || isstring(r.config)
        config_str = char(r.config);

        if isfolder(config_str)
            % Directory: hypothesis dir from decipher output
            fprintf('  Config: hypothesis directory %s\n', config_str);
            build_opts = struct();
            if ~isempty(r.experimental_data)
                build_opts.experimental_data = r.experimental_data;
            end
            if ~isempty(r.organism)
                build_opts.organism = r.organism;
            end
            if ~isempty(r.condition)
                build_opts.condition = r.condition;
            end
            build_opts.output_dir = r.output_dir;
            cfg = build_hypothesis_config(config_str, build_opts);

        elseif isfile(config_str)
            % File path: .m config file
            fprintf('  Config: file %s\n', config_str);
            [cfg_dir, cfg_name, ~] = fileparts(config_str);
            if ~isempty(cfg_dir)
                addpath(cfg_dir);
            end
            cfg_fn = str2func(cfg_name);
            cfg = cfg_fn();
        else
            error('MRHP:validate:ConfigNotFound', ...
                  'Config ''%s'' is not a valid function, file, or directory.', config_str);
        end
    else
        error('MRHP:validate:InvalidConfig', ...
              'Config must be a function handle, file path, or directory path.');
    end

    % Override output_dir
    cfg.output_dir = r.output_dir;

    % Override organism/condition if specified
    if ~isempty(r.organism)
        cfg.organism = r.organism;
    end

    % Overlay experimental data if provided and not already loaded
    if ~isempty(r.experimental_data) && ~isfield(cfg, 'experimental_data')
        if isfolder(r.experimental_data)
            fprintf('  Experimental data dir: %s\n', r.experimental_data);
        end
    end

    fprintf('  Organism:  %s\n', cfg.organism);
    fprintf('  Output:    %s\n', r.output_dir);
    fprintf('  Pipeline:  %s\n\n', r.pipeline_mode);

    % --- Run pipeline ---
    run_pipeline_generic(cfg, r.pipeline_mode);

    fprintf('\n╔══════════════════════════════════════════════════════════╗\n');
    fprintf('║  VALIDATION COMPLETE                                   ║\n');
    fprintf('║  Output: %s\n', r.output_dir);
    fprintf('╚══════════════════════════════════════════════════════════╝\n\n');
end
