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