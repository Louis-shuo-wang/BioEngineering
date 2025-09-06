import numpy as np
import matplotlib.pyplot as plt

# Parameters
D_n = 0.0004608  # diffusion coefficient
chi0 = 0.0599      # chemotactic sensitivity
alpha = 0.6      # saturation constant
dt = 0.00005       # time step
dx = 0.005        # grid spacing

# Domain
N = 50
x = np.linspace(0, 3, N)
y = np.linspace(0, 3, N)
X, Y = np.meshgrid(x, y)

# Synthetic concentration field: Gaussian centered at (1.5, 1.5)
c = np.exp(-0.05 * ((X - 1.5)**2 + (Y - 1.5)**2))
n = np.ones_like(c)  # uniform endothelial density for demonstration

# Compute half-step fluxes F_x and F_y
F_x = np.zeros_like(c)
F_y = np.zeros_like(c)

for j in range(1, N-1):
    for i in range(1, N-1):
        # Interpolated chi and n
        chi = chi0/(1 + alpha * c[j, i])
        n_xp = n[j, i]
        n_yp = n[j, i]

        F_x[j, i] = chi * n[j, i] * (c[j, i+1] - c[j, i]) / dx
        F_y[j, i] = chi * n[j, i] * (c[j+1, i] - c[j, i]) / dx

# Compute directional probabilities on interior (avoid edges)
c_center = c[1:-1, 1:-1]
c_right  = c[1:-1, 2:]
c_left   = c[1:-1, :-2]
c_up     = c[:-2, 1:-1]
c_down   = c[2:, 1:-1]

term_right = (1 / (1 + alpha * c_center) + 1 / (1 + alpha * c_right))
term_left  = (1 / (1 + alpha * c_center) + 1 / (1 + alpha * c_left))
term_up    = (1 / (1 + alpha * c_center) + 1 / (1 + alpha * c_up))
term_down  = (1 / (1 + alpha * c_center) + 1 / (1 + alpha * c_down))

# Adjust Xh, Yh for plotting (interior mesh)
Xh, Yh = X[1:-1, 1:-1], Y[1:-1, 1:-1]

# Plotting
fig1, ax1 = plt.subplots()
for xi in range(4):
    ax1.plot([xi, xi], [0, 3], color='gray', linewidth=0.8)
for yi in range(4):
    ax1.plot([0, 3], [yi, yi], color='gray', linewidth=0.8)
labels = {'P0': (1.5, 1.5), 'P2': (2.5, 1.5), 'P1': (0.5, 1.5),
          'P4': (1.5, 2.5), 'P3': (1.5, 0.5)}
for lbl, (px, py) in labels.items():
    ax1.scatter(px, py, s=100, color='red' if lbl == 'P0' else 'blue')
    ax1.text(px, py - 0.1, lbl, ha='center', va='top', fontsize=15)
ax1.set_aspect('equal')
ax1.axis('off')
fig1.suptitle(
    'Motility Probability Assignment\non a von Neumann Neighborhood',
    fontsize=16
)

# === Plot 2: Chemotactic Flux Field ===
step = 3

fig2, ax2 = plt.subplots()
ax2.quiver(
    X[::step, ::step], Y[::step, ::step],
    F_x[::step, ::step], F_y[::step, ::step],
    angles='xy', scale_units='xy', scale=0.5, width=0.005, color='blue'
)
ax2.set_xlim(0, 3)
ax2.set_ylim(0, 3)
ax2.set_aspect('equal')
ax2.set_title(
    r'Directional Validation of the Chemotactic Flux Field $\boldsymbol{J}_{\mathrm{chemo}}$',
    fontsize=16
)
ax2.tick_params(axis='both', labelsize=15)

plt.tight_layout()
plt.show()

