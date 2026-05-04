"""Quantum-dot / single-electron-transistor energy-level diagram.

Replaces the original green rectangular electrodes with soft-edged red glowing
ellipses (matching the reference shape), and recolors the central dot to black.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle


# ---------- helpers ---------------------------------------------------------
def soft_ellipse(ax, cx, cy, w, h, color=(0.80, 0.05, 0.05),
                 n_rings=60, peak_alpha=0.9, gamma=2.0):
    """Soft-edged red-glow ellipse, built by stacking translucent ellipses.

    The outermost rings fade to zero; the core is near-opaque.  `gamma`
    controls how fast the falloff is (higher → sharper core, softer edge).
    """
    from matplotlib.patches import Ellipse
    for i in range(n_rings, 0, -1):
        frac  = i / n_rings
        alpha = peak_alpha * (1 - frac) ** gamma / n_rings * 6  # cumulative effect
        alpha = max(0.0, min(alpha, 0.25))
        ax.add_patch(Ellipse((cx, cy), w * frac, h * frac,
                             facecolor=color, edgecolor='none', alpha=alpha))
    # Bright core
    ax.add_patch(Ellipse((cx, cy), w * 0.12, h * 0.12,
                         facecolor=(1.0, 0.25, 0.25), edgecolor='none', alpha=0.9))


def tunnel_arrow(ax, x0, y0, x1, y1, rad=-0.5, color='k'):
    """Curved arrow arching above the straight line from (x0,y0) to (x1,y1)."""
    arr = FancyArrowPatch((x0, y0), (x1, y1),
                          connectionstyle=f"arc3,rad={rad}",
                          arrowstyle='-|>', mutation_scale=14,
                          lw=1.3, color=color, zorder=5)
    ax.add_patch(arr)


# ---------- build the figure -------------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 3.8))

# Pale-green panel background with a green border (matches reference)
bg_color     = '#EAF4EA'
border_color = '#6FA36F'
ax.set_facecolor(bg_color)
for spine in ax.spines.values():
    spine.set_edgecolor(border_color)
    spine.set_linewidth(1.4)

# ---- electrodes: soft horizontal red glows (matching the reference shape) ---
left_cx,  right_cx = 1.65, 8.35
lead_cy            = 3.00
lead_w, lead_h     = 3.2, 1.25
soft_ellipse(ax, left_cx,  lead_cy, lead_w, lead_h)
soft_ellipse(ax, right_cx, lead_cy, lead_w, lead_h)

# ---- central quantum-dot energy levels ---
# μ_{N+1}  (the occupied level, with the electron as a black dot)
mu_level_y   = 3.00
mu_above_y   = 4.05
mu_x0, mu_x1 = 4.55, 5.45

ax.plot([mu_x0, mu_x1], [mu_level_y, mu_level_y], color='k', lw=1.8, zorder=4)
ax.plot([mu_x0, mu_x1], [mu_above_y, mu_above_y], color='k', lw=1.8, zorder=4)

# Black electron dot on the lower level
ax.plot((mu_x0 + mu_x1) / 2, mu_level_y, 'o',
        color='k', markersize=10, markeredgecolor='k', zorder=6)

# Level labels (boxed so they read clearly against the background)
bbox = dict(facecolor=bg_color, edgecolor='none', pad=1.0)
ax.text(mu_x1 + 0.12, mu_level_y, r'$\mu_{N+1}$',
        fontsize=14, va='center', ha='left', bbox=bbox, zorder=7)
ax.text((mu_x0 + mu_x1)/2, mu_above_y + 0.32, r'$\mu_{N+2}$',
        fontsize=14, va='center', ha='center', bbox=bbox, zorder=7)

# ---- tunnelling arrows (curved, arching above the barriers) ---
# Left electrode → μ_{N+2} on the dot
tunnel_arrow(ax, left_cx + 0.35, lead_cy + 0.75,
                 mu_x0 - 0.05,  mu_above_y,       rad=-0.45)
# μ_{N+2} → right electrode
tunnel_arrow(ax, mu_x1 + 0.05,  mu_above_y,
                 right_cx - 0.35, lead_cy + 0.75, rad=-0.45)

# ---- V_G gate symbol (capacitor under the dot) ---
gate_x = (mu_x0 + mu_x1) / 2
ax.plot([gate_x, gate_x], [mu_level_y - 0.25, 1.55], color='k', lw=1.3)
# Two capacitor plates
ax.plot([gate_x - 0.30, gate_x + 0.30], [1.40, 1.40], color='k', lw=1.8)
ax.plot([gate_x - 0.22, gate_x + 0.22], [1.20, 1.20], color='k', lw=1.8)
# Gate lead going down
ax.plot([gate_x, gate_x], [1.20, 0.80], color='k', lw=1.3)
# V_G label
ax.text(gate_x + 0.40, 1.30, r'$V_G$', fontsize=14, va='center')

# ---- frame ---
ax.set_xlim(0, 10)
ax.set_ylim(0, 6)
ax.set_aspect('equal')
ax.set_xticks([]); ax.set_yticks([])

plt.tight_layout()
plt.show()
