% This function is used to produce the bifurcation diagram
clear all;clf;

params = initializeParams;

Dn = params.Dn;
Dc = params.Dc;
alpha = params.alpha;
xic = params.xic;
chi0 = params.chi0;
n0 = params.n0;

La = 1;
Lb = 1;

eta_ns = logspace(-10,0, 2001);
m = 20;
wave_vecs = [];
for i = 0:20
    for j = 0:i
        wave_vecs = [wave_vecs ; i,j];
    end
end

max_res = zeros(length(eta_ns),1);
for i = 1:length(eta_ns)
    eta_n = eta_ns(i);
    chi_eff = chi0 / (1 + alpha * eta_n * n0 / xic);
    lambda = (wave_vecs(:, 1)*pi/La).^2 + (wave_vecs(:, 2)*pi/Lb).^2;
    a1 = (Dn + Dc)*lambda + xic;
    a0 = Dn * Dc * lambda.^2 + (Dn * xic - eta_n * chi_eff * n0) * lambda;
    max_re = -inf;
    for j = 1:length(lambda)
        poly = [1, a1(j), a0(j)];
        r = roots(poly);
        max_re = max(max_re, max(real(r)));
    end
    max_res(i) = max_re;
end

set(groot, 'defaultAxesFontSize', 18);
set(groot, 'defaultTextFontSize', 18);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure(1);

loglog(eta_ns, max_res, 'b-');
xlabel('$\eta_n$ (TAF production rate)', 'Interpreter', 'latex');
ylabel('$\max_{\lambda} \Re(\sigma_{+}(\eta_n))$', 'Interpreter', 'latex');
title('Bifurcation Diagram', 'Interpreter', 'latex', 'FontSize', 22);
