function [y_pred, H] = bridge_from_phi(phi, t, lam, beta, y_max)
%BRIDGE_FROM_PHI Compute bridge: y(t) = y_max * (1 - exp(-H(t)^beta))
%   H(t) = integral(lam * phi(t) dt)
%   [y_pred, H] = BRIDGE_FROM_PHI(phi, t, lam, beta, y_max)

if nargin < 5, y_max = 100.0; end

integrand = lam * phi(:)';
H = cumtrapz(t(:)', integrand);
H = max(H, 0);
y_pred = y_max * (1 - exp(-H.^beta));
y_pred = min(max(y_pred, 0), y_max);
end
