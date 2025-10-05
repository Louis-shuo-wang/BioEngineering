clear; close all;
mat1 = load('stripe.mat');
nsnaps1 = mat1.nsnapshots;
csnaps1 = mat1.csnapshots;

La = 4;
Lb = 1;

dx = 0.02;
Nx = La / dx;
x = linspace(0, La, Nx);
y = linspace(0, Lb, Nx);

t_save = [10,20,30];
n_save = numel(t_save);

figure('Position',[100 100 1200 700]); 
tiledlayout(3, 2, 'TileSpacing','tight','Padding','tight');

for k = 1:2*n_save
    nexttile;
    if mod(k,2) == 1
    imagesc(x, y, nsnaps1{1,(k+1)/2}');
    axis equal tight;
    colorbar;
    letter = char('a'+k-1);
    title(sprintf('(%s) n(t=%d)', letter, t_save((k+1)/2)));
    else
       imagesc(x, y, csnaps1{1,k/2}');
       axis equal tight;
       colorbar;
       letter = char('a'+k-1);
       title(sprintf('(%s) c(t=%d)', letter, t_save(k/2)));
    end
end

saveas(gcf, 'Figure 9', 'svg');

mat2 = load('finite_no.mat');
nsnaps2 = mat2.nsnapshots;
csnaps2 = mat2.csnapshots;

mat3 = load('finite_bidirection.mat');
nsnaps3 = mat3.nsnapshots;
csnaps3 = mat3.csnapshots;

mat4 = load('finite_unidirection.mat');
nsnaps4 = mat4.nsnapshots;
csnaps4 = mat4.csnapshots;

mat5 = load('finite_suppression.mat');
nsnaps5 = mat5.nsnapshots;
csnaps5 = mat5.csnapshots;

La = 1;
Lb = 1;

dx = 0.02;
Nx = La / dx;
x = linspace(0, La, Nx);
y = linspace(0, Lb, Nx);

t_save = [10, 20, 30];
n_save = numel(t_save);

figure('Position', [100, 200, 1050, 1800]);
tiledlayout(4,3, 'TileSpacing', 'tight', 'Padding', 'tight');

for k = 1:12
    nexttile;
    if k <= 3
    imagesc(x, y, nsnaps2{1,k}');
    axis equal tight;
    colorbar;
    letter = char('a'+k-1);
    title(sprintf('(%s) n(t=%d)', letter, t_save(k)));
    elseif k <= 6
        imagesc(x, y, nsnaps3{1,k-3}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) n(t=%d)', letter, t_save(k-3)));
    elseif k <= 9
        imagesc(x, y, nsnaps4{1,k-6}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) n(t=%d)', letter, t_save(k-6)));
    else
        imagesc(x, y, nsnaps5{1,k-9}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) n(t=%d)', letter, t_save(k-9)));
    end
end

saveas(gcf, 'Figure 10', 'svg');

figure('Position', [100, 200, 1050, 1800]);
tiledlayout(4,3, 'TileSpacing', 'tight', 'Padding', 'tight');

for k = 1:12
    nexttile;
    if k <= 3
    imagesc(x, y, csnaps2{1,k}');
    axis equal tight;
    colorbar;
    letter = char('a'+k-1);
    title(sprintf('(%s) c(t=%d)', letter, t_save(k)));
    elseif k <= 6
        imagesc(x, y, csnaps3{1,k-3}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) c(t=%d)', letter, t_save(k-3)));
    elseif k <= 9
        imagesc(x, y, csnaps4{1,k-6}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) c(t=%d)', letter, t_save(k-6)));
    else
        imagesc(x, y, csnaps5{1,k-9}');
        axis equal tight;
        colorbar;
        letter = char('a'+k-1);
        title(sprintf('(%s) c(t=%d)', letter, t_save(k-9)));
    end
end

saveas(gcf, 'Figure 11', 'svg');