function gem = load_gem_tsv(gem_dir)
%LOAD_GEM_TSV Load a GEM from three TSV files exported from COBRA/COBRApy.
%
%   gem = LOAD_GEM_TSV(gem_dir)
%
%   Reads a genome-scale metabolic model from three tab-separated files:
%
%     gem_dir/
%       reactions.tsv     — reaction_id, name, equation, gene_rule, ec,
%                           subsystem, reversible, lb, ub
%       metabolites.tsv   — metabolite_id, name, formula, compartment, charge
%       S_matrix.tsv      — rows=metabolites, cols=reactions, integer coefficients
%                           first column header: metabolite\reaction
%
%   Output struct fields:
%     gem.S              — sparse (n_mets × n_rxns) stoichiometric matrix
%     gem.rxn_ids        — cell array of reaction IDs
%     gem.rxn_names      — cell array of reaction names
%     gem.rxn_equations  — cell array of reaction equation strings
%     gem.rxn_genes      — cell array of gene rule strings
%     gem.rxn_ec         — cell array of EC numbers
%     gem.rxn_subsystems — cell array of subsystem annotations
%     gem.rxn_reversible — logical vector (1=reversible)
%     gem.rxn_lb         — lower bounds (n_rxns × 1)
%     gem.rxn_ub         — upper bounds (n_rxns × 1)
%     gem.met_ids        — cell array of metabolite IDs
%     gem.met_names      — cell array of metabolite names
%     gem.met_formulas   — cell array of molecular formulas
%     gem.met_compartments — cell array of compartment tags
%     gem.met_charges    — numeric vector of charges
%     gem.n_rxns         — number of reactions
%     gem.n_mets         — number of metabolites
%
%   Example:
%     gem = load_gem_tsv('inputs/shewanella/gem_tsv/');
%     fprintf('Loaded %d reactions, %d metabolites\n', gem.n_rxns, gem.n_mets);
%
%   See also: DECIPHER_ROUTES, VALIDATE_ROUTE_GAPS

    % --- Validate input directory ---
    if ~isfolder(gem_dir)
        error('MRHP:load_gem_tsv:DirNotFound', ...
              'GEM directory not found: %s', gem_dir);
    end

    rxn_file = fullfile(gem_dir, 'reactions.tsv');
    met_file = fullfile(gem_dir, 'metabolites.tsv');
    smat_file = fullfile(gem_dir, 'S_matrix.tsv');

    for f = {rxn_file, met_file, smat_file}
        if ~isfile(f{1})
            error('MRHP:load_gem_tsv:FileNotFound', ...
                  'Required file not found: %s', f{1});
        end
    end

    % --- Load reactions.tsv ---
    fprintf('  Loading reactions from %s ...\n', rxn_file);
    rxn_raw = readtable(rxn_file, 'FileType','text', 'Delimiter','\t', ...
                        'TextType','string', 'VariableNamingRule','preserve');
    gem.rxn_ids        = cellstr(rxn_raw{:,1});
    gem.rxn_names      = cellstr(rxn_raw{:,2});
    gem.rxn_equations  = cellstr(rxn_raw{:,3});
    gem.rxn_genes      = cellstr(rxn_raw{:,4});
    gem.rxn_ec         = cellstr(rxn_raw{:,5});
    if width(rxn_raw) >= 6
        gem.rxn_subsystems = cellstr(rxn_raw{:,6});
    else
        gem.rxn_subsystems = repmat({''}, size(gem.rxn_ids));
    end
    if width(rxn_raw) >= 7
        gem.rxn_reversible = logical(rxn_raw{:,7});
    else
        gem.rxn_reversible = false(numel(gem.rxn_ids), 1);
    end
    if width(rxn_raw) >= 9
        gem.rxn_lb = rxn_raw{:,8};
        gem.rxn_ub = rxn_raw{:,9};
    else
        gem.rxn_lb = zeros(numel(gem.rxn_ids), 1);
        gem.rxn_ub = 1000 * ones(numel(gem.rxn_ids), 1);
    end
    gem.n_rxns = numel(gem.rxn_ids);

    % --- Load metabolites.tsv ---
    fprintf('  Loading metabolites from %s ...\n', met_file);
    met_raw = readtable(met_file, 'FileType','text', 'Delimiter','\t', ...
                        'TextType','string', 'VariableNamingRule','preserve');
    gem.met_ids          = cellstr(met_raw{:,1});
    gem.met_names        = cellstr(met_raw{:,2});
    gem.met_formulas     = cellstr(met_raw{:,3});
    gem.met_compartments = cellstr(met_raw{:,4});
    if width(met_raw) >= 5
        gem.met_charges = met_raw{:,5};
    else
        gem.met_charges = zeros(numel(gem.met_ids), 1);
    end
    gem.n_mets = numel(gem.met_ids);

    % --- Load S_matrix.tsv ---
    fprintf('  Loading stoichiometric matrix from %s ...\n', smat_file);
    smat_raw = readtable(smat_file, 'FileType','text', 'Delimiter','\t', ...
                         'TextType','string', 'VariableNamingRule','preserve');

    % First column = metabolite IDs (row labels)
    smat_met_ids = cellstr(smat_raw{:,1});
    % Remaining columns = reaction coefficients
    smat_rxn_ids = smat_raw.Properties.VariableNames(2:end);
    S_dense = smat_raw{:, 2:end};
    if isa(S_dense, 'string')
        S_dense = double(S_dense);
    end
    S_dense = double(S_dense);
    gem.S = sparse(S_dense);

    % --- Validate consistency ---
    % S_matrix rows must match metabolites.tsv
    if size(gem.S, 1) ~= gem.n_mets
        warning('MRHP:load_gem_tsv:MetMismatch', ...
                'S_matrix has %d rows but metabolites.tsv has %d entries. Using S_matrix row labels.', ...
                size(gem.S, 1), gem.n_mets);
        gem.met_ids = smat_met_ids;
        gem.n_mets = numel(gem.met_ids);
    end

    % S_matrix columns must match reactions.tsv
    if size(gem.S, 2) ~= gem.n_rxns
        warning('MRHP:load_gem_tsv:RxnMismatch', ...
                'S_matrix has %d columns but reactions.tsv has %d entries. Using S_matrix column labels.', ...
                size(gem.S, 2), gem.n_rxns);
        gem.rxn_ids = smat_rxn_ids(:);
        gem.n_rxns = numel(gem.rxn_ids);
    end

    fprintf('  GEM loaded: %d metabolites × %d reactions (%d non-zero entries)\n', ...
            gem.n_mets, gem.n_rxns, nnz(gem.S));
end
