import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

df = pd.read_csv('../mass_loss_CO10.csv')

stars      = df['Source'].tolist()
Mdot_thin  = df['Mdot_thin_Msunyr'].astype(float).values
Mdot_thick = df['Mdot_thick_Msunyr'].astype(float).values
T_mb       = df['T_mb_peak_K'].astype(float).values
types      = df['Type'].tolist()

color_map = {
    'C star': '#378ADD',
    'S star': '#1D9E75',
    'Mira M': '#D85A30',
}

T_norm = (T_mb - T_mb.min()) / (T_mb.max() - T_mb.min() + 1e-9)
sizes  = 80 + T_norm * 320
colors = [color_map.get(t, '#888780') for t in types]
x      = np.arange(len(stars))


fig, ax = plt.subplots(figsize=(10, 6))

ax.scatter(x, Mdot_thin,  s=sizes, c=colors, marker='o', alpha=0.85, zorder=3)
ax.scatter(x, Mdot_thick, s=sizes, c=colors, marker='s', alpha=0.5,  zorder=3)

for i in range(len(stars)):
    ax.plot([x[i], x[i]], [Mdot_thin[i], Mdot_thick[i]],
            color=colors[i], linewidth=1, alpha=0.4, zorder=2)

ax.axhline(2e-6, color='gray', linestyle='--', linewidth=1,
           alpha=0.6, label='Lit. RV Aqr (2×10⁻⁶)')


ax.set_yscale('log')
ax.set_xticks(x)
ax.set_xticklabels(stars, fontsize=11, rotation=20, ha='right')
ax.set_ylabel(r'$\dot{M}$ (M$_\odot$ yr$^{-1}$)', fontsize=13)
ax.set_title('Mass loss rate CO(1-0)', fontsize=14, fontweight='bold')
ax.grid(True, which='both', alpha=0.2, linestyle='--')
ax.tick_params(axis='both', direction='in', which='both',
               length=5, width=1.2, top=True, right=True)
for spine in ax.spines.values():
    spine.set_linewidth(1.3)

type_patches = [mpatches.Patch(color=c, label=t)
                for t, c in color_map.items()
                if t in types]  # solo los tipos que aparecen en el CSV

# marcadores thin vs thick
marker_thin  = plt.scatter([], [], marker='o', color='gray', s=80,
                           label='Ṁ opticamente delgado')
marker_thick = plt.scatter([], [], marker='s', color='gray', s=80,
                           alpha=0.5, label='Ṁ opticamente grueso')

# tamaño por T_mb — 3 valores representativos del rango real
t_vals = np.linspace(T_mb.min(), T_mb.max(), 3)
size_patches = [
    plt.scatter([], [], marker='o', color='gray',
                s=80 + ((t - T_mb.min()) / (T_mb.max() - T_mb.min() + 1e-9)) * 320,
                label=f'T_mb = {t:.1f} K')
    for t in t_vals
]

legend1 = ax.legend(handles=type_patches + size_patches,
                    title='Tipo / T_mb', fontsize=9,
                    loc='upper right', framealpha=0.9)
ax.add_artist(legend1)
ax.legend(handles=[marker_thin, marker_thick],
          fontsize=9, loc='lower right', framealpha=0.9)

plt.tight_layout()
plt.savefig('mass_loss_CO10.png', dpi=150, bbox_inches='tight')
plt.show()