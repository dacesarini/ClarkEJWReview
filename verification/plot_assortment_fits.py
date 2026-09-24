"""Create the observed-versus-fitted PNAS education comparison figure."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

HERE = Path(__file__).resolve().parent
report = json.loads((HERE / 'assortment_fit_results.json').read_text(encoding='utf-8'))
fits = [f for f in report['fits'] if f['cohort'] == '1860-1919' and f['weighting'] == 'unweighted']
genetic = next(f for f in fits if f['model'] == 'genetic')
ph_fixed = next(f for f in fits if f['model'] == 'phenotypic' and not f['free_theta'])
ph_free = next(f for f in fits if f['model'] == 'phenotypic' and f['free_theta'])
order = [1, 0, 3, 2, 4, 5, 6, 7, 8, 9, 10]
labels = ['Father–son', 'Full brothers', 'Grandfather–grandson', 'Uncle–nephew',
          'First cousins', 'First cousins, once removed', 'Second cousins',
          'Second cousins, once removed', 'Third cousins',
          'Third cousins, once removed', 'Fourth cousins']
observed = np.array(genetic['observed'])[order]
gpred = np.array(genetic['predictions'])[order]
blue, orange, ink = '#2369A0', '#C16B25', '#202E3B'
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11,
                     'axes.labelcolor': ink, 'text.color': ink,
                     'xtick.color': '#586775', 'ytick.color': ink,
                     'pdf.fonttype': 42, 'ps.fonttype': 42})
fig, axes = plt.subplots(1, 2, figsize=(13.8, 8.6), sharey=True)
fig.subplots_adjust(left=.235, right=.975, top=.735, bottom=.19, wspace=.13)
fig.patch.set_facecolor('white')
fig.text(.035, .954, 'Which model matches the observed family resemblance?',
         fontsize=21, weight='bold', va='top')
fig.text(.035, .910, 'Higher education · English male birth cohort, 1860–1919 · PNAS Table 2',
         fontsize=12.5, color='#586775', va='top')
legend = [Line2D([], [], color=ink, marker='D', linestyle='None', markersize=7, label='Observed correlation'),
          Line2D([], [], color=blue, marker='o', linestyle='None', markersize=8, label='Phenotypic assortment'),
          Line2D([], [], color=orange, marker='s', linestyle='None', markersize=7, label='Genetic assortment')]
fig.legend(handles=legend, loc='upper left', bbox_to_anchor=(.03,.868),
           frameon=False, ncol=3, columnspacing=2.4, handletextpad=.5)

for ax, fit, title in zip(axes, [ph_fixed, ph_free], ['θ = 1  |  No measurement error', 'θ free  |  Measurement error allowed']):
    ppred = np.array(fit['predictions'])[order]
    for row in (0, 2):
        ax.axhspan(row-.48, row+.48, facecolor='#F0F4F7', zorder=0)
    for row, obs, p, g in zip(range(11), observed, ppred, gpred):
        ax.hlines(row, min(obs,p,g), max(obs,p,g), color='#B8C2CB', linewidth=1.1, zorder=2)
        ax.vlines(obs, row-.19, row+.19, color='#CAD1D7', linewidth=.8, zorder=2)
    yy = np.arange(11)
    ax.scatter(ppred, yy-.15, s=60, c=blue, edgecolors='white', linewidths=.6, zorder=4)
    ax.scatter(gpred, yy+.15, s=51, c=orange, marker='s', edgecolors='white', linewidths=.6, zorder=4)
    ax.scatter(observed, yy, s=38, c=ink, marker='D', edgecolors='white', linewidths=.5, zorder=5)
    ax.set_title(title, fontsize=13, weight='bold', loc='left', pad=34)
    ax.text(0, 1.028, f"Squared error: phenotypic {fit['raw_sse']:.4f}  ·  genetic {genetic['raw_sse']:.4f}",
            transform=ax.transAxes, fontsize=10, color='#586775')
    ax.set_xlim(0,.48)
    ax.set_ylim(10.65,-.65)
    ax.set_xticks(np.arange(0,.5,.1))
    ax.set_xlabel('Correlation', labelpad=10)
    ax.set_yticks(yy, labels)
    ax.tick_params(axis='both', length=0, pad=9)
    ax.grid(axis='x', color='#E6EBEF', linewidth=.8, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

fig.text(.035,.109, 'Shaded rows are lineal relatives. Small vertical offsets separate the markers; only horizontal position measures correlation.',
         fontsize=10, color='#586775')
fig.text(.035,.077, 'Fits minimize unweighted squared errors in correlation levels across all 11 relationships. No observed spouse moment is included.',
         fontsize=10, color='#586775')
fig.text(.035,.045, 'Genetic predictions are unchanged when θ is freed: only θ × h² is identified. Observed correlations are estimates, not known population truths.',
         fontsize=10, color='#586775')
destination = HERE.parent / 'figures'
destination.mkdir(exist_ok=True)
for extension in ('png','pdf','svg'):
    path = destination / f'assortment_fits_education_1860_1919.{extension}'
    fig.savefig(path, dpi=200, facecolor='white')
    print(path)
plt.close(fig)
