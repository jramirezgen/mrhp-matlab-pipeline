function [lam, beta_val, y_max, R2, RMSE, y_pred] = fit_bridge(phi, t, D_ref, y_max_init)
%FIT_BRIDGE Fit lambda, beta, y_max to match D_ref via global optimization.
%   [lam, beta_val, y_max, R2, RMSE, y_pred] = FIT_BRIDGE(phi, t, D_ref, y_max_init)

if nargin < 4 || isempty(y_max_init)
    y_max_init = max(D_ref);
    if y_max_init <= 0, y_max_init = 95.0; end
end

lb = [0.001, 0.1, max(1e-6, y_max_init*0.3)];
ub = [500.0,  5.0, y_max_init*2.0];

% Objective
    function sse = obj(params)
        y_p = bridge_from_phi(phi, t, params(1), params(2), params(3));
        sse = sum((D_ref(:)' - y_p).^2);
    end

% Multi-start optimization
best_sse = Inf;
best_x = [1, 1, y_max_init];
rng(42);
n_starts = 50;
for i = 1:n_starts
    x0_trial = lb + rand(1,3) .* (ub - lb);
    try
        opts_opt = optimset('Display','off','MaxFunEvals',5000,'MaxIter',2000,'TolFun',1e-10);
        [x_opt, f_opt] = fminsearch(@obj, x0_trial, opts_opt);
        % Project into bounds
        x_opt = max(min(x_opt, ub), lb);
        f_opt = obj(x_opt);
        if f_opt < best_sse
            best_sse = f_opt;
            best_x = x_opt;
        end
    catch
    end
end

lam = best_x(1);
beta_val = best_x(2);
y_max = best_x(3);

y_pred = bridge_from_phi(phi, t, lam, beta_val, y_max);
ss_res = sum((D_ref(:)' - y_pred).^2);
ss_tot = sum((D_ref(:)' - mean(D_ref)).^2);
R2 = 1 - ss_res / max(ss_tot, 1e-12);
RMSE = sqrt(mean((D_ref(:)' - y_pred).^2));
end
