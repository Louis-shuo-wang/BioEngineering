clear all; close all;

params = initializeParams;

Dn = 100 * params.Dn;
Dc = 0.5 * params.Dc;
alpha = params.alpha;
xic = 5 * params.xic;
chi0 = 20 * params.chi0;
n0 = params.n0;
eta_n = 1;

La = 1;
Lb = 1;

dx = 0.02;
Nx = La / dx;
x = linspace(0, La, Nx);
y = linspace(0, Lb, Nx);

dt = 5e-3;
Tfinal = 100;
Nt = ceil(Tfinal/dt);
t_save = [10,20,30];
n_save = numel(t_save);

nsnapshots = cell(1,n_save);
csnapshots = cell(1,n_save);


Tfinal = t_save(end);
Nt = Tfinal/dt;
time = (1:Nt)'*dt;

p = eta_n;
c0 = eta_n * n0 /xic;

rng(1);
[X,Y] = meshgrid(x,y);
n_mat = n0 * (1 + 1e-3 * randn(Nx,Nx));
c_mat = c0 * (1 + 1e-3 * randn(Nx,Nx));

n = n_mat(:);
c = c_mat(:);

n_snap = zeros(Nx,Nx,n_save);
c_snap = zeros(Nx,Nx,n_save);
save_idx = 1;

%% 2D Laplacian

L2D = laplacian2D(Nx, dx);
I_N = speye(Nx*Nx);

M_n_L = I_N - (dt/2)*Dn*L2D;
M_n_R = I_N + (dt/2)*Dn*L2D;

M_c_L = I_N - (dt/2)*Dc*L2D;
M_c_R = I_N + (dt/2)*Dc*L2D;

% diagnostics: compute critical p and k2
chi_eff_fun = @(pval) chi0 ./ (1 + alpha * pval / xic);
pc_fun = @(pval) pval * chi_eff_fun(pval) * n0 - Dn * xic;

chi_eff_val = chi_eff_fun(p);
num = p * chi_eff_val * n0 - Dn * xic;
num1 = num/(Dn*Dc);
if num > 0
    k2 = sqrt(num / (Dn * Dc));
    k_star = pi/k2;
else
    k2 = NaN; k_star = NaN;
end
fprintf('Diagnostics: p=%.3g, chi_eff=%.3g, Dn*xi_c=%.3g, p*chi_eff*n0=%.3g\n', ...
    p, chi_eff_val, Dn*xic, p*chi_eff_val*n0);
fprintf('k2 = %.4g, wave = %.4f, critical k* = %.4f', k2, num1,k_star);


% diagnostics: compute critical p and 
for step = 1:Nt
    t = step*dt;
    n_mat = reshape(n,Nx,Nx);
    c_mat = reshape(c,Nx,Nx);

 c_pad_x = [c_mat(:,2),c_mat,c_mat(:,end-1)];
    c_pad_y = [c_mat(2,:);c_mat;c_mat(end-1,:)];
    c_grad_x = (c_pad_x(:,3:end)-c_pad_x(:,1:end-2))/(2*dx);
    c_grad_y = (c_pad_y(3:end,:)-c_pad_y(1:end-2,:))/(2*dx);

    % --- chemotactic velocities
    chi = chi0 ./ (1+alpha*c_mat);
    u_x = chi .* c_grad_x;
    u_y = chi .* c_grad_y;

    % --- neighbors of n (for upwind)
    nL = [n_mat(:,1), n_mat(:,1:end-1)];   % left
    nR = [n_mat(:,2:end), n_mat(:,end)];   % right
    nD = [n_mat(1,:); n_mat(1:end-1,:)];   % down
    nU = [n_mat(2:end,:); n_mat(end,:)];   % up

    % --- upwind flux in x
    F_nx = zeros(size(n_mat));
    pos = (u_x > 0);
    F_nx(pos)  = u_x(pos) .* nL(pos);
    F_nx(~pos) = u_x(~pos) .* nR(~pos);

    % --- upwind flux in y
    F_ny = zeros(size(n_mat));
    pos = (u_y > 0);
    F_ny(pos)  = u_y(pos) .* nD(pos);
    F_ny(~pos) = u_y(~pos) .* nU(~pos);

    % --- divergence of fluxes
    F_nx_pad = [F_nx(:,2), F_nx, F_nx(:,end-1)];
    F_ny_pad = [F_ny(2,:); F_ny; F_ny(end-1,:)];
    div_n = (F_nx_pad(:,3:end) - F_nx_pad(:,1:end-2))/(2*dx) + ...
            (F_ny_pad(3:end,:) - F_ny_pad(1:end-2,:))/(2*dx);

    fn = -div_n(:);

    fc = (p * n_mat - xic * c_mat);
    fc = fc(:);

    n = M_n_L \ (M_n_R * n + dt * fn);
    c = M_c_L \ (M_c_R * c + dt * fc);

    n = max(n,0);
    c = max(c,0);

    if save_idx <= n_save && abs(t-t_save(save_idx))<dt/2
        n_snap(:,:,save_idx) = reshape(n,Nx,Nx);
        c_snap(:,:,save_idx) = reshape(c,Nx,Nx);
        nsnapshots{1,save_idx} = n_snap(:,:,save_idx);
        csnapshots{1,save_idx} = c_snap(:,:,save_idx);
        fprintf('Saved snapshot at t=%.3f (step %d)\n', t, step);
        save_idx = save_idx + 1;
    end
end

set(groot, 'defaultAxesFontSize', 12);
set(groot, 'defaultTextFontSize', 12);
set(groot, 'defaultAxesLineWidth', 1.2);
set(groot, 'defaultLineLineWidth', 2);
set(groot, 'defaultAxesTickDir', 'out');
set(groot, 'defaultAxesBox', 'off');

figure('Position',[100 100 1800 500]); 
tiledlayout(1, n_save, 'TileSpacing','compact','Padding','compact');

for k = 1:n_save
    nexttile;
    imagesc(x, y, n_snap(:,:,k)');
    axis equal tight;
    title(sprintf('n, t=%.1f', t_save(k)));
end
colorbar;
sgtitle('Snapshots of n (Neumann BC)');

figure('Position',[100 100 1800 500]); 
tiledlayout(1, n_save, 'TileSpacing','compact','Padding','compact');

for k = 1:n_save
    nexttile;
    imagesc(x, y, c_snap(:,:,k)');
    axis equal tight;
    title(sprintf('c, t=%.1f', t_save(k)));
end
colorbar;
sgtitle('Snapshots of c (Neumann BC)');

save("finite_suppression.mat","params","nsnapshots","csnapshots")