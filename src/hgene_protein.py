"""Shared AD169 protein normalization for database preparation and reporting.

HGVS event/3' rules, three-letter keys and the project's one-letter p. display.
This is deliberately not a general HGVS parser or a cross-strain liftover tool.
"""
from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

CMV_ACCESSION = "FJ527563.1"
CMV_FASTA_SHA256 = "2ffdb56b0aa949e473aab53941e0e124c9477e5b12b81f9a8d8e6f6fe6da41c8"
AA = "ACDEFGHIKLMNPQRSTVWY*"
AA3 = dict(zip(AA, (
    "Ala", "Cys", "Asp", "Glu", "Phe", "Gly", "His", "Ile", "Lys", "Leu",
    "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr", "Ter",
)))
HGVS_VERSION = "1"
HGVS_COLUMNS = ("hgvs_p", "hgvs_status", "hgvs_note")

# Code génétique standard ; * indique un codon stop.
CODONS = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

@dataclass(frozen=True)
class Normalization:
    raw: str
    normalized: str = ""
    status: str = "UNSUPPORTED"
    note: str = ""
    kind: str = ""

    @property
    def ok(self) -> bool:
        return self.status == "OK"

    @property
    def hgvs_p(self) -> str:
        """Unparenthesized, three-letter comparison key; never a partial event."""
        return re.sub(r"[A-Z*]", lambda m: AA3[m[0]], self.normalized) if self.ok else ""


class ProteinError(ValueError):
    def __init__(self, status: str, message: str):
        self.status = status
        super().__init__(message)


def translate(sequence: str) -> str:
    return "".join(CODONS.get(sequence[i:i + 3], "X")
                   for i in range(0, len(sequence) - 2, 3))


def aa_offset(chrom: str) -> int:
    # UL89_3 is the first 888 nt, UL89_1 the last 1137 nt in HHV5.fasta.
    return 296 if chrom == "UL89_1-HHV5" else 0


def protein_gene(value: object) -> str:
    text = str(value or "").strip().upper()
    match = re.search(r"\b(UL\d+[A-Z]?)", text)
    return match[1] if match else re.sub(r"[^A-Z0-9]+", "", text)


class CMVReference:
    def __init__(self, fasta: Path):
        self.fasta = fasta
        self.sha256 = hashlib.sha256(fasta.read_bytes()).hexdigest()
        if self.sha256 != CMV_FASTA_SHA256:
            raise ValueError("FASTA is not the validated hgene AD169/FJ527563.1 reference")
        self.dna: dict[str, str] = {}
        key = ""
        for line in fasta.read_text().splitlines():
            if line.startswith(">"):
                key = line[1:].split()[0]
                self.dna[key] = ""
            elif line.strip():
                self.dna[key] += line.strip().upper()
        self.proteins = {k.split("-")[0]: translate(v) for k, v in self.dna.items()
                         if not k.startswith("UL89_")}
        self.proteins["UL89"] = translate(
            self.dna["UL89_3-HHV5"] + self.dna["UL89_1-HHV5"])

    def vcf_problem(self, metadata: dict[str, str]) -> str:
        accession = metadata.get("hgene_ref_accession", "")
        digest = metadata.get("hgene_ref_fasta_sha256", "")
        if accession and accession != CMV_ACCESSION:
            return f"VCF reference {accession!r} differs from {CMV_ACCESSION}"
        if digest and digest != self.sha256:
            return "VCF reference FASTA SHA-256 differs from validated AD169"
        if not accession and not digest:
            return "VCF has neither hgene_ref_accession nor hgene_ref_fasta_sha256"
        return ""

    def check_dna(self, chrom: str, pos: int, ref: str) -> bool:
        dna = self.dna.get(chrom, "")
        return pos > 0 and bool(ref) and dna[pos - 1:pos - 1 + len(ref)] == ref


def hgvs_metadata(reference: CMVReference) -> dict[str, str]:
    return {"hgvs_reference": CMV_ACCESSION,
            "hgvs_fasta_sha256": reference.sha256,
            "hgvs_normalizer_version": HGVS_VERSION}


def normalize_database_mutation(record: Mapping[str, str], reference: CMVReference,
                                coordinates: str = "global") -> Normalization:
    """Offline only: current curated BCSQ is authoritative, never *_original."""
    raw = record.get("aa_chg_bcsq") or record.get("aa_change") or ""
    gene = protein_gene(record.get("gene") or record.get("chrom"))
    chrom = record.get("chrom", "")
    coordinates = record.get("aa_coordinate_system") or coordinates
    if chrom and (not chrom.endswith("-HHV5") or protein_gene(chrom) != gene):
        return Normalization(raw, status="GENE_MISMATCH", note="Gene/contig is not from this CMV reference")
    if coordinates not in {"global", "local"} or (coordinates == "local" and gene == "UL89" and not chrom):
        return Normalization(raw, status="AMBIGUOUS", note="UL89 local coordinates require a segment contig")
    if record.get("ref_accession") and record["ref_accession"] != CMV_ACCESSION:
        return Normalization(raw, status="REFERENCE_MISMATCH", note="Database reference differs from AD169")
    offset = aa_offset(chrom) if coordinates == "local" else 0
    return normalize_protein(raw, reference.proteins.get(gene, ""), offset)


def normalize_variant_mutation(chrom: str, gene: str, pos: int, dna_ref: str,
                               raw: str, consequence: str, reference: CMVReference) -> Normalization:
    """Online: normalize one BCSQ consequence against the same full protein."""
    gene = protein_gene(gene)
    if gene != protein_gene(chrom):
        return Normalization(raw, status="GENE_MISMATCH", note="BCSQ gene differs from contig")
    if not reference.check_dna(chrom, pos, dna_ref):
        return Normalization(raw, status="DNA_REF_MISMATCH", note="VCF REF differs from AD169")
    return normalize_protein(raw, reference.proteins.get(gene, ""), aa_offset(chrom), consequence)


def load_prepared_records(records: Sequence[Mapping[str, str]], metadata: Mapping[str, str],
                          reference: CMVReference) -> list[dict[str, str]]:
    """Check provenance and use saved keys, without renormalizing the database."""
    if not records:
        return []
    for key, expected in hgvs_metadata(reference).items():
        if metadata.get(key) != expected:
            raise ValueError(f"Database {key} missing or different; run hgene_prepare_resistance_db.py --fasta first")
    result = []
    for record in records:
        if any(key not in record for key in HGVS_COLUMNS):
            raise ValueError("Database needs hgvs_p, hgvs_status and hgvs_note; run hgene_prepare_resistance_db.py first")
        gene = protein_gene(record.get("gene") or record.get("chrom"))
        status = record.get("hgvs_status", "")
        key = record.get("hgvs_p", "")
        note = record.get("hgvs_note", "")
        if status == "OK" and not key.startswith("p."):
            raise ValueError(f"Empty/invalid hgvs_p at database row {record.get('source_row', '?')}; rerun preprocessing")
        if not status:
            raise ValueError("Database contains an unprepared row; rerun hgene_prepare_resistance_db.py")
        chrom = record.get("chrom", "")
        if gene not in reference.proteins or (chrom and (protein_gene(chrom) != gene or not chrom.endswith("-HHV5"))):
            status, note = "GENE_MISMATCH", "Database gene/contig differs from the prepared CMV reference"
        result.append(dict(record, normalization_gene=gene, normalization_status=status,
                           normalization_note=note, normalization_reference=CMV_ACCESSION,
                           normalization_fasta_sha256=reference.sha256))
    return result


def normalize_protein(raw: str, protein: str, offset: int = 0,
                      consequence: str = "") -> Normalization:
    """Normalize one definite event; ambiguous inputs retain raw text and a status."""
    original = raw
    raw = raw.strip().removeprefix("p.")
    if raw.startswith("(") and raw.endswith(")"):
        raw = raw[1:-1]
    aa1 = {three: one for one, three in AA3.items()}
    raw = re.sub("|".join(aa1), lambda m: aa1[m[0]], raw)
    if not protein:
        return Normalization(original, status="NO_REFERENCE", note="Unknown gene")

    def check(pos: int, residues: str) -> None:
        if pos < 1 or pos + len(residues) - 1 > len(protein):
            raise ProteinError("OUT_OF_RANGE", f"Position {pos} outside reference")
        if protein[pos - 1:pos - 1 + len(residues)] != residues:
            raise ProteinError("REF_MISMATCH", f"Reference residues differ at {pos}")

    def residue(pos: int) -> str:
        if not 1 <= pos <= len(protein):
            raise ProteinError("OUT_OF_RANGE", f"Position {pos} outside reference")
        return protein[pos - 1]

    try:
        synonymous = re.fullmatch(r"(\d+)([A-Z*])", raw)
        if synonymous and "synonymous" in consequence:
            check(int(synonymous.group(1)) + offset, synonymous.group(2))
            return Normalization(original, "p.=", "OK", kind="synonymous")
        # Legacy single-AA deletion, e.g. 981D>981del.
        deletion = re.fullmatch(r"(\d+)([A-Z*])>(\d+)del", raw)
        if deletion:
            p, aa, q = deletion.groups()
            if p != q:
                raise ProteinError("AMBIGUOUS", "Different deletion coordinates")
            raw = f"{aa}{p}del"

        # BCSQ describes a reference and alternate block at the same start.
        block = re.fullmatch(r"(\d+)([A-Z*]+)>(?:(\d+))?([A-Z*]+)", raw)
        if block:
            p, removed, q, inserted = block.groups()
            if q and int(p) != int(q):
                raise ProteinError("AMBIGUOUS", "Different BCSQ block starts")
            start = int(p) + offset - 1
            check(start + 1, removed)
            if "frameshift" in consequence:
                # A BCSQ block need not provide the entire new reading frame.
                i = 0
                while i < min(len(removed), len(inserted)) and removed[i] == inserted[i]:
                    i += 1
                if i >= min(len(removed), len(inserted)):
                    raise ProteinError("INCOMPLETE", "First altered frameshift residue unavailable")
                pos = start + i + 1
                if inserted[i] == "*":
                    return Normalization(original, f"p.{removed[i]}{pos}*", "OK", kind="stop_gained")
                stop = inserted.find("*", i)
                suffix = f"*{stop - i + 1}" if stop >= 0 else ""
                return Normalization(original, f"p.{removed[i]}{pos}{inserted[i]}fs{suffix}",
                                     "PARTIAL", "Frameshift: no protein-equivalence matching", "frameshift")
            end = start + len(removed)
        else:
            fs = re.fullmatch(r"([A-Z])(\d+)([A-Z])?fs(?:\*(\d+|\?))?", raw)
            if fs:
                aa, p, alt, stop = fs.groups()
                pos = int(p) + offset
                check(pos, aa)
                return Normalization(original, f"p.{aa}{pos}{alt or ''}fs" + (f"*{stop}" if stop else ""),
                                     "PARTIAL", "Frameshift: no protein-equivalence matching", "frameshift")
            sub = re.fullmatch(r"([A-Z*])(\d+)([A-Z*]|=)", raw)
            event = re.fullmatch(r"([A-Z*]?)(\d+)(?:_([A-Z*]?)(\d+))?(delins|del|dup|ins)([A-Z*]*)", raw)
            if sub:
                aa, p, inserted = sub.groups()
                start = int(p) + offset - 1
                check(start + 1, aa)
                end = start + 1
                if inserted == "=":
                    inserted = aa
            elif event:
                left, p, right, q, kind, inserted = event.groups()
                p = int(p) + offset
                q = int(q) + offset if q else p
                if p > q:
                    raise ProteinError("AMBIGUOUS", "Reversed coordinate range")
                check(p, left or residue(p))
                check(q, right or residue(q))
                if kind == "ins":
                    if q != p + 1 or not inserted:
                        raise ProteinError("AMBIGUOUS", "Insertion needs two adjacent flanking positions")
                    start = end = p
                elif kind == "dup":
                    if inserted:
                        raise ProteinError("UNSUPPORTED", "Duplication with extra sequence")
                    inserted = protein[p - 1:q]
                    start = end = q
                else:
                    if (kind == "del" and inserted) or (kind == "delins" and not inserted):
                        raise ProteinError("UNSUPPORTED", "Malformed deletion/delins")
                    start, end = p - 1, q
            else:
                status = "AMBIGUOUS" if "in>" in raw or "ins" in raw else "UNSUPPORTED"
                raise ProteinError(status, "Notation requires manual review")

        if set(inserted) - set(AA):
            raise ProteinError("UNSUPPORTED", "Unknown amino acid")
        if "*" in protein[start:end]:
            raise ProteinError("UNSUPPORTED", "Stop loss/extension requires additional sequence")
        if start == 0:
            raise ProteinError("UNSUPPORTED", "Initiation-codon change requires manual review")
        removed = protein[start:end]
        if set(removed) - set(AA):
            raise ProteinError("UNSUPPORTED", "Unknown reference amino acid")
        if "*" in inserted and len(removed) == 1 and inserted == "*":
            return Normalization(original, f"p.{removed}{start + 1}*", "OK", kind="stop_gained")
        if "*" in inserted:
            raise ProteinError("UNSUPPORTED", "Complex termination event requires manual review")
        alternate = protein[:start] + inserted + protein[end:]
        if alternate == protein:
            return Normalization(original, "p.=", "OK", kind="synonymous")

        # Longest common prefix FIRST gives the rightmost equivalent minimal
        # event, including rotated tandem motifs. Then strip the common suffix.
        left = 0
        while left < min(len(protein), len(alternate)) and protein[left] == alternate[left]:
            left += 1
        right = 0
        while (right < min(len(protein) - left, len(alternate) - left)
               and protein[len(protein) - right - 1] == alternate[len(alternate) - right - 1]):
            right += 1
        end = len(protein) - right
        inserted = alternate[left:len(alternate) - right if right else len(alternate)]
        removed = protein[left:end]

        def label(a: int, b: int) -> str:
            return f"{residue(a)}{a}" + (f"_{residue(b)}{b}" if a != b else "")

        if not removed:
            size = len(inserted)
            if left >= size and protein[left - size:left] == inserted:
                text, kind = f"{label(left - size + 1, left)}dup", "duplication"
            elif 0 < left < len(protein):
                text, kind = f"{label(left, left + 1)}ins{inserted}", "insertion"
            else:
                raise ProteinError("UNSUPPORTED", "Terminal insertion")
        elif not inserted:
            text, kind = f"{label(left + 1, end)}del", "deletion"
        elif len(removed) == len(inserted) == 1:
            text, kind = f"{label(left + 1, end)}{inserted}", "substitution"
        else:
            text, kind = f"{label(left + 1, end)}delins{inserted}", "delins"
        assert protein[:left] + inserted + protein[end:] == alternate
        return Normalization(original, "p." + text, "OK", kind=kind)
    except ProteinError as exc:
        return Normalization(original, status=exc.status, note=str(exc))


def main(argv=None) -> int:
    """Small TSV bridge: the R report uses the same Python normalizer as Word."""
    import argparse
    import csv
    import sys

    parser = argparse.ArgumentParser(description="Normalize CMV report variants and check prepared HGVS database keys")
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--variants", required=True, type=Path, help="TSV with CHROM, POS, REF, gene, amino_acid_change, consequence")
    parser.add_argument("--database", required=True, type=Path, help="prepared database as TSV")
    parser.add_argument("--metadata", required=True, type=Path, help="database metadata TSV with key/value columns")
    parser.add_argument("--output-variants", required=True, type=Path)
    parser.add_argument("--output-database", required=True, type=Path)
    parser.add_argument("--vcf-reference", default="")
    parser.add_argument("--vcf-fasta-sha256", default="")
    args = parser.parse_args(argv)

    def read_table(path):
        with path.open(encoding="utf-8-sig", newline="") as handle:
            return [{k.lower(): (v or "") for k, v in row.items()} for row in csv.DictReader(handle, delimiter="\t")]

    def write_table(path, rows, columns):
        with path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    try:
        protected = {p.resolve() for p in (args.fasta, args.variants, args.database, args.metadata)}
        outputs = {args.output_variants.resolve(), args.output_database.resolve()}
        if len(outputs) != 2 or outputs & protected:
            raise ValueError("Output tables must be distinct from each other and from inputs")
        reference = CMVReference(args.fasta)
        problem = reference.vcf_problem({"hgene_ref_accession": args.vcf_reference,
                                         "hgene_ref_fasta_sha256": args.vcf_fasta_sha256})
        if problem:
            raise ValueError(problem)
        metadata = {r["key"]: r["value"] for r in read_table(args.metadata)}
        records = load_prepared_records(read_table(args.database), metadata, reference)
        db_columns = ["normalization_gene", "normalization_status", "normalization_note"]
        db_rows = [{k: r[k] for k in db_columns} for r in records]
        rows = []
        for v in read_table(args.variants):
            norm = normalize_variant_mutation(v["chrom"], v["gene"], int(v["pos"]), v["ref"],
                                              v["amino_acid_change"], v["consequence"], reference)
            rows.append(dict(hgvs_p=norm.hgvs_p, normalization_gene=protein_gene(v["gene"]),
                             normalization_status=norm.status, normalization_note=norm.note,
                             aa_normalized=norm.normalized))
        write_table(args.output_variants, rows, ["hgvs_p", "normalization_gene", "normalization_status",
                                                "normalization_note", "aa_normalized"])
        write_table(args.output_database, db_rows, db_columns)
    except (OSError, ValueError, KeyError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())