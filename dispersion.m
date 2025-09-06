% This function is used to produce dispersion diagram.

clear all; clf;
params = initializeParams;

Dn = params.Dn;
Dc = params.Dc;
alpha = params.alpha;
xic = params.xic;
chi0 = params.chi0;
n0 = params.n0;

La = 1;
Lb = 1;

eta_ns = [0.001, 0.01, 0.05, 0.1, 0.5, 1];
lambda_Upper = 0;
for i = 1:length(eta_ns)
    eta_n = eta_ns(i);
    chi_eff = chi0 / (1+alpha * eta_n * n0 / xic);
    lambda_upper = (eta_n * chi_eff * n0 - Dn * xic) / (Dn * Dc);
    lambda_Upper = max(lambda_upper, lambda_Upper);
end

colors = lines(length(eta_ns)+1);
linestyles = {'-', '--', '-', '--', '-', '--'};

max_res = zeros(length(eta_ns), 101);
for i = 1:length(eta_ns)
    eta_n = eta_ns(i);
    chi_eff = chi0 / (1+alpha * eta_n * n0 / xic);
lambda = linspace(0, lambda_Upper, 101);
a0 = Dn * Dc *  lambda .^ 2 + (Dn * xic - eta_n * chi_eff * n0) * lambda;
a1 = (Dn + Dc) * lambda + xic;
max_re = [];
for j = 1:length(lambda)
    p = [1 a1(j) a0(j)];
    r = roots(p);
    max_re = [max_re, max(real(r))];
end
max_res(i,:) = max_re;
end

set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure(1);
for i = 1:length(eta_ns)
    label_str = sprintf('$\\eta_n = %.1e$', eta_ns(i));
    plot(lambda, max_res(i, :), 'Color', colors(i,:), 'LineStyle', linestyles{i}, 'DisplayName', label_str); hold on;
end

yline(0, 'Color', colors(end,:), 'LineWidth', 2, 'HandleVisibility', 'off');

xLimits = xlim;
yLimits = ylim;

text(0.5 * xLimits(2), 0.03 * (yLimits(2) - yLimits(1)), '$\Re(\sigma)=0$', ...
    'Interpreter','latex', 'Color', colors(end,:), ...
    'HorizontalAlignment', 'right');
xlabel('$\lambda$', 'Interpreter', 'latex');
ylabel('$\Re(\sigma_{+}(\lambda))$', 'Interpreter', 'latex');
legend("show", 'Interpreter', 'latex');
title('Dispersion curves', 'Interpreter', 'latex');
