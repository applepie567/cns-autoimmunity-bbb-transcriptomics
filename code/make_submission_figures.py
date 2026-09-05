#!/usr/bin/env python3
"""Regenerate all final figures from the frozen derived tables."""

import subprocess
import sys
from pathlib import Path


def main() -> None:
    code = Path(__file__).resolve().parent
    for script in [
        "make_bbi_figures.py",
        "make_orthogonal_supplementary_figure.py",
        "make_s9_from_frozen_results.py",
        "build_supplementary_figures_pdf.py",
    ]:
        subprocess.run([sys.executable, str(code / script)], check=True)


if __name__ == "__main__":
    main()
