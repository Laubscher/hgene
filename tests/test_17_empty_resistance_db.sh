#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

PYTHONPATH="$REPO_ROOT/src${PYTHONPATH:+:$PYTHONPATH}" python3 - "$REPO_ROOT" <<'PY'
import subprocess
import sys
import tempfile
from pathlib import Path

from docx import Document
from openpyxl import Workbook

root = Path(sys.argv[1])
db = root / "data/db/hsv1/resistances.xlsx"

from hgene_clinical_report import (
    EMPTY_RESISTANCE_DB_REQUIRED_COLUMNS,
    build_conclusion,
    read_resistance_db,
)

records, metadata, readable = read_resistance_db(db)
assert readable, "the complete empty template should be readable"
assert records == [], "the default database should contain no resistance records"
assert "db_name" in metadata
assert build_conclusion([], [], {}, True, bool(records), True) == "", \
    "an empty database must not produce a susceptible conclusion"

with tempfile.TemporaryDirectory(prefix="hgene-empty-db-test-") as temp:
    temp = Path(temp)
    template = temp / "template.docx"
    doc = Document()
    doc.add_paragraph("_MUT-ACV-UL23_")
    doc.add_paragraph("_INT-ACV-UL23_")
    doc.save(template)
    vcf = temp / "empty.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
        encoding="utf-8",
    )
    coverage = temp / "coverage.tsv"
    coverage.write_text("gene\tlow_coverage_20x\tlow_coverage_100x\n", encoding="utf-8")
    report = temp / "report.docx"
    variants_tsv = temp / "clinical_variants.tsv"
    run = subprocess.run(
        [sys.executable, str(root / "src/hgene_clinical_report.py"),
         "--template", str(template), "--output", str(report),
         "--clinical-variants-tsv", str(variants_tsv),
         "--vcf", str(vcf), "--coverage", str(coverage),
         "--resistance-db", str(db), "--virus", "HHV1"],
        text=True,
        capture_output=True,
        check=False,
    )
    assert run.returncode == 0, run.stderr
    assert "resistance database is empty" in run.stderr, run.stderr
    assert "clinical variants TSV requires" not in run.stderr, run.stderr
    assert report.is_file(), "clinical report should still be written"
    assert variants_tsv.read_text(encoding="utf-8").startswith("gene\tchrom\tpos\t")

with tempfile.TemporaryDirectory(prefix="hgene-invalid-empty-db-test-") as temp:
    incomplete = Path(temp) / "incomplete.xlsx"
    wb = Workbook()
    sheet = wb.active
    sheet.title = "resistances"
    sheet.append(sorted(EMPTY_RESISTANCE_DB_REQUIRED_COLUMNS - {"hgvs_note"}))
    sheet.append(["NA"] * sheet.max_column)
    wb.save(incomplete)
    _, _, readable = read_resistance_db(incomplete)
    assert not readable, "an empty workbook missing a required field is not valid"

print("OK test_17_empty_resistance_db")
PY
