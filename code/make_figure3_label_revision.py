"""Replot frozen Figure 3 with readable panel A labels, without reanalysis."""
from pathlib import Path
import importlib.util
import json
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.text import Text

R = Path(__file__).resolve().parents[1]
OUT = R / 'output' / 'figures'
SOURCE = R / 'baseline' / 'BBI_v2.0.0' / 'code' / 'make_bbi_figures.py'


def revised_save(fig, stem):
    # The three small axes following the panel A title are the donor profiles.
    labels = ['Control', 'Chronic active', 'Chronic inactive']
    profile_axes = fig.axes[1:4]
    for ax in profile_axes:
        assert len(ax.collections) == 3
        ax.set_xticks(range(3), labels, rotation=90, ha='center', va='top')
        ax.tick_params(axis='x', labelsize=6.8, pad=3)
    # Match the manuscript's black lettering, including numbers on dark cells.
    for item in fig.findobj(match=Text):
        if item.get_color() == 'white':
            item.set_path_effects([pe.withStroke(linewidth=.65, foreground='white')])
        item.set_color('black')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    ticks = [t for ax in profile_axes for t in ax.get_xticklabels()]
    boxes = [t.get_window_extent(renderer) for t in ticks]
    assert all(not a.overlaps(b) for i,a in enumerate(boxes) for b in boxes[i+1:])
    next_title = fig.axes[5]._left_title.get_window_extent(renderer)
    # Reserve 12 points above the following panel heading for the longer labels.
    clearance = 12 * fig.dpi / 72
    needed = next_title.y1 + clearance - min(b.y0 for b in boxes)
    if needed > 0:
        shift = needed / fig.bbox.height
        for ax in profile_axes:
            pos = ax.get_position()
            ax.set_position([pos.x0, pos.y0 + shift, pos.width, pos.height - shift])
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        boxes = [t.get_window_extent(renderer) for t in ticks]
    assert min(b.y0 for b in boxes) > next_title.y1 + 3
    OUT.mkdir(exist_ok=True, parents=True)
    name = 'Figure_3_human_MS_microenvironment_stage_revised'
    for ext in ['png', 'pdf']:
        target = OUT / (name + '.' + ext)
        temp = OUT / (name + '.tmp.' + ext)
        fig.savefig(temp, dpi=350, bbox_inches='tight', facecolor='white', format=ext)
        os.replace(temp, target)
    plt.close(fig)
    print(json.dumps({'figure': 3, 'panel': 'A', 'tick_labels': labels,
                      'rotation_degrees': 90, 'label_overlap': False,
                      'source_results': 'unchanged'}))


if __name__ == '__main__':
    spec = importlib.util.spec_from_file_location('frozen_bbi_figures', SOURCE)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.set_style()
    matplotlib.rcParams['pdf.fonttype'] = 42
    module.save = revised_save
    module.fig3_human_ms()
