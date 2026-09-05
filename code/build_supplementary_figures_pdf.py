#!/usr/bin/env python3
"""Assemble the current Supplementary Figures S1–S9 into one PDF."""

from pathlib import Path

from PIL import Image
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas


ROOT = Path(__file__).resolve().parents[1]
FIGURE_DIR = ROOT / "figures" / "supplementary"
OUTPUT = ROOT / "figures" / "BBI_Supplementary_Figures_S1-S9.pdf"

FIGURES = [
    ("S1", "Figure_S1_GSE210776_QC.png"),
    ("S2", "Figure_S2_acute_EAE_sample_scores.png"),
    ("S3", "Figure_S3_GSE199460_subtype_sensitivity.png"),
    ("S4", "Figure_S4_GSE95401_time_course.png"),
    ("S5", "Figure_S5_GSE279183_EC_QC.png"),
    ("S6", "Figure_S6_GSE208747_threshold_sensitivity.png"),
    ("S7", "Figure_S7_VEGFA_all_modules.png"),
    ("S8", "Figure_S8_published_functional_protein_source_data.png"),
    ("S9", "Figure_S9_GSE284005_spatial_validation.png"),
]


def main() -> None:
    missing = [name for _, name in FIGURES if not (FIGURE_DIR / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing supplementary figures: " + ", ".join(missing))

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    page_width, page_height = landscape(A4)
    pdf = canvas.Canvas(str(OUTPUT), pagesize=(page_width, page_height), pageCompression=1)
    pdf.setTitle("BBI Supplementary Figures S1–S9")
    pdf.setAuthor("Manuscript supplementary material")

    left = right = 24
    bottom = 20
    header_top = 20
    header_height = 24
    available_width = page_width - left - right
    available_height = page_height - bottom - header_top - header_height

    for page_number, (figure_id, filename) in enumerate(FIGURES, start=1):
        path = FIGURE_DIR / filename
        with Image.open(path) as img:
            pixel_width, pixel_height = img.size

        pdf.bookmarkPage(figure_id)
        pdf.addOutlineEntry(f"Supplementary Figure {figure_id}", figure_id, level=0, closed=False)
        pdf.setFillColorRGB(0.16, 0.19, 0.22)
        pdf.setFont("Helvetica-Bold", 11)
        pdf.drawString(left, page_height - header_top, f"Supplementary Figure {figure_id}")
        pdf.setFont("Helvetica", 8)
        pdf.setFillColorRGB(0.42, 0.46, 0.49)
        pdf.drawRightString(page_width - right, page_height - header_top, f"{page_number} / {len(FIGURES)}")

        scale = min(available_width / pixel_width, available_height / pixel_height)
        draw_width = pixel_width * scale
        draw_height = pixel_height * scale
        x = (page_width - draw_width) / 2
        y = bottom + (available_height - draw_height) / 2
        pdf.drawImage(
            ImageReader(str(path)), x, y, width=draw_width, height=draw_height,
            preserveAspectRatio=True, mask="auto",
        )
        pdf.showPage()

    pdf.save()
    print(OUTPUT)


if __name__ == "__main__":
    main()
