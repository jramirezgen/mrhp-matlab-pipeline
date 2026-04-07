function sto = tau_leaping_expression(u_t, u_vals, ktx_nM, n_cells, base_seed, expr_params)
%TAU_LEAPING_EXPRESSION Stochastic tau-leaping for gene expression.
%   sto = TAU_LEAPING_EXPRESSION(u_t, u_vals, ktx_nM, n_cells, base_seed, expr_params)
%   Returns struct with mean_m, std_m, mean_p, std_p, t.
%
%   Vectorized across cells for performance (vs scalar-per-cell baseline).

if nargin < 4 || isempty(n_cells),    n_cells   = 20;   end
if nargin < 5 || isempty(base_seed),  base_seed = 42;   end
if nargin < 6 || isempty(expr_params), expr_params = struct(); end

tau     = getfield_default(expr_params, 'stochastic_tau', 0.01);
beta_m  = getfield_default(expr_params, 'beta_m', 8.32);
beta_p  = getfield_default(expr_params, 'beta_p', 0.924);
ktl     = getfield_default(expr_params, 'ktl', 6.0);
nm2mol  = getfield_default(expr_params, 'nm_to_mol', 6022);
t_span  = getfield_default(expr_params, 't_span', [0 72]);
t_eval  = getfield_default(expr_params, 't_eval', (0:0.5:72)');

ktx_mol = ktx_nM * nm2mol;

% Time grid
t_max = t_span(2);
n_steps = round(t_max / tau);
t_all = (0:n_steps)' * tau;

% Interpolate u signal (compute once)
u_interp = interp1(u_t(:), u_vals(:), t_all, 'linear', u_vals(end));

% Output grid
t_out = t_eval(:)';
n_out = numel(t_out);

% Precompute output-step mapping
out_step = round(t_out / tau);

% Set single reproducible seed
rng(base_seed);

% State vectors for all cells simultaneously (vectorized)
m_all = zeros(n_cells, 1);
p_all = zeros(n_cells, 1);

% Preallocate output storage
rec_m = zeros(n_cells, n_out);
rec_p = zeros(n_cells, n_out);
out_idx = 1;

for step = 1:n_steps
    u_now = max(u_interp(step), 0);

    % Transcription: Poisson -- vectorized across all cells
    rate_txn = ktx_mol * u_now * tau;
    if rate_txn > 0
        n_txn = poissrnd(repmat(rate_txn, n_cells, 1));
    else
        n_txn = zeros(n_cells, 1);
    end

    % Translation: Poisson -- vectorized (rate depends on per-cell m)
    rate_tln = ktl * m_all * tau;
    n_tln = poissrnd(max(rate_tln, 0));

    % mRNA degradation: Poisson, capped at current count
    rate_deg_m = beta_m * m_all * tau;
    n_deg_m = min(m_all, poissrnd(max(rate_deg_m, 0)));

    % Protein degradation: Poisson, capped
    rate_deg_p = beta_p * p_all * tau;
    n_deg_p = min(p_all, poissrnd(max(rate_deg_p, 0)));

    % Update state
    m_all = max(m_all + n_txn - n_deg_m, 0);
    p_all = max(p_all + n_tln - n_deg_p, 0);

    % Record at output times
    while out_idx <= n_out && out_step(out_idx) <= step
        rec_m(:, out_idx) = m_all;
        rec_p(:, out_idx) = p_all;
        out_idx = out_idx + 1;
    end
end

% Convert molecules -> nM
rec_m_nM = rec_m / nm2mol;
rec_p_nM = rec_p / nm2mol;

sto.t      = t_out(:)';
sto.mean_m = mean(rec_m_nM, 1);
sto.std_m  = std(rec_m_nM, 0, 1);
sto.mean_p = mean(rec_p_nM, 1);
sto.std_p  = std(rec_p_nM, 0, 1);
sto.n_cells = n_cells;
sto.all_m  = rec_m;       % molecules/cell, n_cells x n_out
sto.all_m_nM = rec_m_nM;  % nM, n_cells x n_out
end

function val = getfield_default(s, fname, default)
    if isfield(s, fname), val = s.(fname); else, val = default; end
end
