function [is_connected, gaps, report] = validate_route_gaps(S_sub, species, reactions, target_idx, seed_idx)
%VALIDATE_ROUTE_GAPS Check connectivity of a metabolic sub-network.
%
%   [is_connected, gaps, report] = VALIDATE_ROUTE_GAPS(S_sub, species, ...
%       reactions, target_idx, seed_idx)
%
%   Verifies that a continuous path exists from seed metabolites to the
%   target in the stoichiometric sub-network.  Common cofactors (NAD, NADH,
%   ATP, ADP, CoA, H2O, H+, CO2, Pi) are excluded from gap analysis since
%   they participate in many reactions without being true pathway
%   intermediates.
%
%   Inputs:
%     S_sub       — (n_sp × n_rxn) stoichiometric matrix of the sub-network
%     species     — cell array of species IDs (n_sp × 1)
%     reactions   — cell array of reaction IDs (n_rxn × 1)
%     target_idx  — index of the target metabolite in species
%     seed_idx    — index (or vector of indices) of seed metabolites
%
%   Outputs:
%     is_connected — true if every non-cofactor species is reachable
%     gaps         — struct array with fields:
%                      .metabolite  — species ID of the orphan
%                      .idx         — index in species
%                      .type        — 'produced_only' | 'consumed_only' | 'isolated'
%     report       — cell array of diagnostic strings
%
%   Example:
%     [ok, g, r] = validate_route_gaps(S, species, rxns, 6, 1);
%     if ~ok, fprintf('Gaps found: %d\n', numel(g)); end
%
%   See also: DECIPHER_ROUTES, LOAD_GEM_TSV

    % --- Cofactors to exclude from gap analysis ---
    cofactor_patterns = {'cpd00003', 'cpd00004', 'cpd00006', ... % NAD+, NADH, NADP+
                         'cpd00005', ...                         % NADPH
                         'cpd00002', 'cpd00008', ...             % ATP, ADP
                         'cpd00009', 'cpd00012', ...             % Pi, PPi
                         'cpd00010', ...                         % CoA
                         'cpd00001', ...                         % H2O
                         'cpd00067', ...                         % H+
                         'cpd00011', ...                         % CO2
                         'cpd00013', ...                         % NH3
                         'NAD', 'NADH', 'ATP', 'ADP', 'CoA', ...
                         'H2O', 'H+', 'CO2', 'Pi'};

    n_sp  = size(S_sub, 1);
    n_rxn = size(S_sub, 2);

    % --- Identify cofactors ---
    is_cofactor = false(n_sp, 1);
    for i = 1:n_sp
        sp = species{i};
        for j = 1:numel(cofactor_patterns)
            if contains(sp, cofactor_patterns{j}, 'IgnoreCase', true)
                is_cofactor(i) = true;
                break;
            end
        end
    end

    % --- Build adjacency: metabolite reachability via reactions ---
    % A metabolite is "reachable" if there is a reaction that produces it
    % from another reachable metabolite (or it is a seed).
    reachable = false(n_sp, 1);
    reachable(seed_idx) = true;

    % Iterate until no more species become reachable (BFS-like)
    changed = true;
    max_iter = n_rxn + 1;
    iter = 0;
    while changed && iter < max_iter
        changed = false;
        iter = iter + 1;
        for r = 1:n_rxn
            col = full(S_sub(:, r));
            substrates = find(col < 0);
            products   = find(col > 0);

            % If any substrate is reachable, all products become reachable
            if any(reachable(substrates))
                for p = products(:)'
                    if ~reachable(p)
                        reachable(p) = true;
                        changed = true;
                    end
                end
            end
        end
    end

    % --- Detect gaps ---
    gaps = struct('metabolite', {}, 'idx', {}, 'type', {});
    report = {};
    gap_count = 0;

    for i = 1:n_sp
        if is_cofactor(i), continue; end
        if reachable(i), continue; end

        col_participation = full(S_sub(i, :));
        is_produced = any(col_participation > 0);
        is_consumed = any(col_participation < 0);

        if is_produced && ~is_consumed
            gtype = 'produced_only';
        elseif ~is_produced && is_consumed
            gtype = 'consumed_only';
        else
            gtype = 'isolated';
        end

        gap_count = gap_count + 1;
        gaps(gap_count).metabolite = species{i};
        gaps(gap_count).idx = i;
        gaps(gap_count).type = gtype;
        report{end+1} = sprintf('GAP: %s (idx %d) — %s', species{i}, i, gtype); %#ok<AGROW>
    end

    % --- Check if target is reachable ---
    target_reachable = reachable(target_idx);
    if ~target_reachable
        report{end+1} = sprintf('CRITICAL: target %s not reachable from seeds', species{target_idx});
    end

    is_connected = (gap_count == 0) && target_reachable;

    if is_connected
        report{end+1} = sprintf('PASS: All %d non-cofactor species connected. Target reachable.', ...
                                sum(~is_cofactor));
    else
        report{end+1} = sprintf('FAIL: %d gaps detected. Target reachable: %s', ...
                                gap_count, mat2str(target_reachable));
    end

    fprintf('  Gap validation: %d gaps, target reachable=%d\n', gap_count, target_reachable);
end
