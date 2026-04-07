function [t_out, mRNA, protein] = solve_deterministic_expression(t, u_vals, ktx, beta_m, beta_p, ktl, t_eval)
%SOLVE_DETERMINISTIC_EXPRESSION Deterministic expression ODE.
%   dm/dt = ktx * u(t) - beta_m * m
%   dp/dt = ktl * m    - beta_p * p
%   [t_out, mRNA, protein] = SOLVE_DETERMINISTIC_EXPRESSION(t, u_vals, ktx, beta_m, beta_p, ktl, t_eval)

if nargin < 4 || isempty(beta_m), beta_m = 8.32;  end
if nargin < 5 || isempty(beta_p), beta_p = 0.924;  end
if nargin < 6 || isempty(ktl),    ktl    = 6.0;    end
if nargin < 7 || isempty(t_eval), t_eval = (0:0.5:72)'; end

% Interpolant for u
t = t(:);
u_vals = u_vals(:);
u_interp = griddedInterpolant(t, u_vals, 'linear', 'nearest');

odefun = @(tt, y) [ktx * u_interp(tt) - beta_m * y(1);
                    ktl * y(1)          - beta_p * y(2)];

opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-12, 'NonNegative', [1 2]);
[t_out, y_out] = ode15s(odefun, t_eval, [0; 0], opts);

mRNA    = y_out(:, 1)';   % 1 x n_time
protein = y_out(:, 2)';
t_out   = t_out';
end
