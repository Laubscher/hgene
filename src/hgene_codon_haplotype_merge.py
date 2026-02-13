#!/usr/bin/env python3
"""
hgene_codon_haplotype_merge.py

- For each codon containing 2 or 3 SNV, reconstructs read-level haplotypes (patterns of R/A)
- Emits every non-reference haplotype whose frequency >= hap_af_min (default 0.10)
- AF is recalculated as: HAP_CT / HAP_DP (HAP_DP = informative_reads)
- DP_CODON is min per-site depth across the codon positions (after MAPQ/baseQ filters)
- Adds HAP/HAP_CT/HAP_DP and CODON_REF/CODON_ALT, AA_REF/AA_ALT, MERGED_POS/MERGED_AF, GENE in INFO
- Logs to stderr: haplotype distribution per tested codon; logs when no haplotype >= threshold.

Notes:
- VCF output may contain multiple records for the same codon position (one per haplotype >= threshold).
"""

import sys
import gzip
from collections import defaultdict
from datetime import datetime

try:
    import pysam  # type: ignore
except Exception:
    pysam = None

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]

def opn(path):
    if path == "-" or path is None:
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def load_fasta(path: str) -> dict:
    seqs = {}
    name, buf = None, []
    with opn(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line)
        if name is not None:
            seqs[name] = "".join(buf).upper()
    return seqs

def parse_gff_cds(path: str):
    cds = []
    with opn(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, phase, attrs = cols
            if ftype != "CDS":
                continue
            start = int(start); end = int(end)
            phase = 0 if phase == "." else int(phase)
            gene = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if "=" not in kv:
                    continue
                k, v = kv.split("=", 1)
                if k in ("gene", "Name", "Parent", "ID") and gene is None:
                    gene = v
            cds.append(dict(chrom=chrom, start=start, end=end, strand=strand, phase=phase, gene=gene or "CDS"))
    by_chr = defaultdict(list)
    for c in cds:
        by_chr[c["chrom"]].append(c)
    return by_chr

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

def aa(codon: str) -> str:
    return CODON_TABLE.get(codon.upper().replace("U", "T"), "X")

def parse_info(s: str) -> dict:
    if s in (".", "", None):
        return {}
    d = {}
    for part in s.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            d[k] = v
        else:
            d[part] = True
    return d

def format_info(d: dict) -> str:
    if not d:
        return "."
    out = []
    for k, v in d.items():
        out.append(k if v is True else f"{k}={v}")
    return ";".join(out)

def get_af(info_d, fmt_keys, sample_fields):
    # INFO/AF
    if "AF" in info_d:
        try:
            return float(str(info_d["AF"]).split(",")[0])
        except Exception:
            pass
    # FORMAT/AF
    if fmt_keys and sample_fields and "AF" in fmt_keys:
        try:
            return float(sample_fields[fmt_keys.index("AF")].split(",")[0])
        except Exception:
            pass
    # FORMAT/AD
    if fmt_keys and sample_fields and "AD" in fmt_keys:
        try:
            ad = sample_fields[fmt_keys.index("AD")].split(",")
            if len(ad) >= 2:
                r = float(ad[0]); a = float(ad[1])
                return a / (r + a) if (r + a) > 0 else None
        except Exception:
            pass
    return None

def _read_base_at_refpos(read, ref_pos0: int):
    """Return (base, qpos) at ref_pos0 (0-based), or (None, None) if not covered."""
    qpos_at = None
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos == ref_pos0:
            if qpos is None:
                return None, None
            qpos_at = qpos
            break
    if qpos_at is None:
        return None, None
    base = read.query_sequence[qpos_at]
    # SAM stores SEQ in read orientation; for reverse-strand alignments, complement the base
    if read.is_reverse:
        base = base.translate(DNA_COMP)
    return base, qpos_at

def phase_haplotypes_sites(   #replace by phase_haplotypes_sites_pileup
    bam_path: str,
    chrom: str,
    sites,
    min_mapq: int = 40,
    min_baseq: int = 20,
    min_informative_reads: int = 10,
):
    """
    Reconstruct read-level haplotypes on 2-3 genomic sites (multi-allelic aware).

    sites: list of dicts with keys:
      - pos (1-based genomic)
      - alleles: set of allowed alleles at this pos (REF + all ALTs), A/C/G/T only
      - frame: 0/1/2 within CDS codon (with respect to codon_ref_cds returned by codon_ctx)

    Returns dict:
      - hap_counts: haplotypes as observed base strings like "GA", "GAC", ...
      - informative_reads: reads covering ALL positions with baseQ>=min_baseq (and MAPQ>=min_mapq)
      - dp_pos: per-position valid depth
      - dp_min: min(dp_pos)
      - ok_min_reads: informative_reads >= min_informative_reads
    """
    if pysam is None:
        raise RuntimeError("You need to install pysam.")

    pos0s = [s["pos"] - 1 for s in sites]
    start0 = min(pos0s)
    end0 = max(pos0s)

    bam = pysam.AlignmentFile(bam_path, "rb")
    hap_counts = defaultdict(int)
    dp_pos = [0] * len(sites)
    informative = 0

    for read in bam.fetch(chrom, start0, end0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.query_qualities is None:
            continue

        bases = []
        ok_all_positions = True

        for j, s_ in enumerate(sites):
            base, qpos = _read_base_at_refpos(read, s_["pos"] - 1)
            if base is None or qpos is None:
                ok_all_positions = False
                break
            if read.query_qualities[qpos] < min_baseq:
                ok_all_positions = False
                break

            b = base.upper()
            if b not in "ACGT":
                ok_all_positions = False
                break

            # valid depth at this position for this read
            dp_pos[j] += 1

            if b not in s_["alleles"]:
                ok_all_positions = False
                break

            bases.append(b)

        if not ok_all_positions:
            continue

        informative += 1
        hap_counts["".join(bases)] += 1

    bam.close()
    dp_min = min(dp_pos) if dp_pos else 0

    return {
        "hap_counts": dict(hap_counts),
        "informative_reads": informative,
        "dp_pos": dp_pos,
        "dp_min": dp_min,
        "ok_min_reads": informative >= min_informative_reads,
    }

def phase_haplotypes_sites_pileup(
    bam_path: str,
    chrom: str,
    sites,  # list of dict: {"pos": int(1-based), "alleles": set("ACGT")}
    min_mapq: int = 40,
    min_baseq: int = 20,
    min_informative_reads: int = 10,
):
    import pysam

    # positions 1-based -> set for fast lookup
    pos_set = {s["pos"] for s in sites}
    sites_by_pos = {s["pos"]: s for s in sites}
    start0 = min(pos_set) - 1
    end0 = max(pos_set)     # pysam pileup end is 0-based exclusive typically; we’ll be safe with end0

    bam = pysam.AlignmentFile(bam_path, "rb")

    # read_name -> {pos: base}
    per_read = defaultdict(dict)

    # Build per-read bases at each target pos using pileup (fast C core)
    for col in bam.pileup(
        chrom,
        start0,
        end0,
        truncate=True,
        stepper="samtools",
        min_mapping_quality=min_mapq,
    ):
        pos1 = col.reference_pos + 1
        if pos1 not in pos_set:
            continue

        for pr in col.pileups:
            if pr.is_del or pr.is_refskip:
                continue
            qpos = pr.query_position
            if qpos is None:
                continue

            aln = pr.alignment
            # baseQ
            if aln.query_qualities is None:
                continue
            if aln.query_qualities[qpos] < min_baseq:
                continue

            b = aln.query_sequence[qpos].upper()
            if b not in "ACGT":
                continue

            # only accept bases that are in REF/ALT set for this site
            if b not in sites_by_pos[pos1]["alleles"]:
                continue

            per_read[aln.query_name][pos1] = b

    bam.close()

    # Count haplotypes across sites (ordered by position)
    ordered_pos = sorted(pos_set)
    hap_counts = defaultdict(int)
    informative = 0
    dp_pos = [0] * len(ordered_pos)

    for _, d in per_read.items():
        # must cover all positions
        if len(d) != len(ordered_pos):
            continue
        informative += 1
        hap = "".join(d[p] for p in ordered_pos)
        hap_counts[hap] += 1
        for j, p in enumerate(ordered_pos):
            dp_pos[j] += 1

    dp_min = min(dp_pos) if dp_pos else 0
    return {
        "hap_counts": dict(hap_counts),
        "informative_reads": informative,
        "dp_pos": dp_pos,
        "dp_min": dp_min,
        "ok_min_reads": informative >= min_informative_reads,
        "ordered_pos": ordered_pos,
    }




def codon_ctx(chrom, pos1, fasta, cds_by_chr):
    """
    Returns (cds, codon_anchor, frame, codon_ref_CDS)
    - pos1: 1-based genomic
    - codon_anchor: + => codon start (1-based)
                    - => rightmost base of triplet on genome (1-based)
    """
    for cds in cds_by_chr.get(chrom, []):
        if cds["start"] <= pos1 <= cds["end"]:
            strand = cds["strand"]
            phase = cds["phase"]
            seq = fasta.get(chrom)
            if not seq:
                return None
            if strand == "+":
                idx = (pos1 - cds["start"]) - phase
                if idx < 0:
                    return None
                frame = idx % 3
                codon_start = pos1 - frame
                codon = seq[codon_start-1:codon_start-1+3]
                if len(codon) != 3:
                    return None
                return (cds, codon_start, frame, codon)
            else:
                idx = (cds["end"] - pos1) - phase
                if idx < 0:
                    return None
                frame = idx % 3
                codon_anchor = pos1 + frame
                triplet = seq[codon_anchor-3:codon_anchor]
                if len(triplet) != 3:
                    return None
                codon = revcomp(triplet)
                return (cds, codon_anchor, frame, codon)
    return None

def codon_alt_from_bases(codon_ref_cds: str, strand: str, sites, hap_bases: str) -> str:
    """
    Build CODON_ALT in CDS orientation by applying observed bases at each site.
    - hap_bases: observed bases across sites in genomic orientation (A/C/G/T), same order as sites
    - sites: each has 'frame' (0/1/2) with respect to codon_ref_cds
    """
    cod = list(codon_ref_cds)
    for s, b in zip(sites, hap_bases):
        bb = b.upper()
        if strand == "-":
            bb = revcomp(bb)
        cod[s["frame"]] = bb
    return "".join(cod)

def ensure_info_headers(header_lines):
    def has(id_):
        return any(h.startswith("##INFO=") and f"ID={id_}," in h for h in header_lines)

    defs = [
        ("MERGED_POS", '##INFO=<ID=MERGED_POS,Number=1,Type=String,Description="Original SNV positions merged into one codon variant">'),
        ("MERGED_AF",  '##INFO=<ID=MERGED_AF,Number=1,Type=String,Description="AFs of original SNVs (same order as MERGED_POS)">'),
        ("CODON_REF",  '##INFO=<ID=CODON_REF,Number=1,Type=String,Description="Reference codon (CDS orientation)">'),
        ("CODON_ALT",  '##INFO=<ID=CODON_ALT,Number=1,Type=String,Description="Mutated codon (CDS orientation)">'),
        ("AA_REF",     '##INFO=<ID=AA_REF,Number=1,Type=String,Description="Reference amino acid">'),
        ("AA_ALT",     '##INFO=<ID=AA_ALT,Number=1,Type=String,Description="Mutated amino acid">'),
        ("GENE",       '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene/Parent/ID extracted from GFF">'),
        ("DP_CODON",   '##INFO=<ID=DP_CODON,Number=1,Type=Integer,Description="Minimum per-site depth across merged codon positions (after MAPQ/baseQ filters)">'),
        ("HAP",        '##INFO=<ID=HAP,Number=1,Type=String,Description="Observed bases across merged positions (same order as MERGED_POS)">'),
        ("HAP_CT",     '##INFO=<ID=HAP_CT,Number=1,Type=Integer,Description="Read count supporting this haplotype pattern">'),
        ("HAP_DP",     '##INFO=<ID=HAP_DP,Number=1,Type=Integer,Description="Informative read depth used to estimate haplotype AF">'),
    ]
    extras = [line for (id_, line) in defs if not has(id_)]

    out = []
    inserted = False
    for h in header_lines:
        if (not inserted) and h.startswith("#CHROM"):
            out.extend(extras)
            inserted = True
        out.append(h)
    if not inserted:
        out.extend(extras)
    return out

def consequence_from_aa(aa_ref: str, aa_alt: str) -> str:
    if aa_ref == aa_alt:
        return "synonymous"
    if aa_alt == "*":
        return "stop_gained"
    if aa_ref == "*":
        return "stop_lost"
    return "missense"

def update_bcsq_value(bcsq: str, new_consequence: str, aa_pos: int,
                      aa_ref: str, aa_alt: str, nuc_change: str) -> str:
    """
    Update first BCSQ record:
      consequence|gene|...|protein_coding|+|88W>88R|262TGG>CGT
    Replace:
      [0] consequence
      [5] AA field
      [6] nuc field
    """
    if not bcsq or bcsq == ".":
        return bcsq
    entries = bcsq.split(",")
    p = entries[0].split("|")
    if len(p) < 7:
        return bcsq
    p[0] = new_consequence
    p[5] = f"{aa_pos}{aa_ref}>{aa_pos}{aa_alt}"
    p[6] = nuc_change
    entries[0] = "|".join(p)
    return ",".join(entries)

def main(vcf_path, fasta_path, gff_path, bam_path,
         hap_af_min=0.10, min_informative_reads=10, min_mapq=40, min_baseq=20):
    fasta = load_fasta(fasta_path)
    cds_by_chr = parse_gff_cds(gff_path)

    header = []
    recs = []

    with opn(vcf_path) as f:
        for line in f:
            if line.startswith("#"):
                header.append(line.rstrip("\n"))
                continue
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4].split(",")[0]
            info_d = parse_info(fields[7] if len(fields) > 7 else ".")
            fmt_keys = fields[8].split(":") if len(fields) > 8 and fields[8] != "." else None
            sample_fields = fields[9].split(":") if len(fields) > 9 and fmt_keys else None
            af = get_af(info_d, fmt_keys, sample_fields)

            ctx = None
            if len(ref) == 1 and len(alt) == 1 and chrom in fasta:
                ctx = codon_ctx(chrom, pos, fasta, cds_by_chr)

            recs.append({
                "line": line,
                "fields": fields,
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "info_d": info_d,
                "af": af,
                "ctx": ctx,  # (cds, codon_anchor, frame, codon_ref_cds)
            })

    # Group SNVs by codon
    groups = defaultdict(list)
    for i, r in enumerate(recs):
        if r["ctx"] is None:
            continue
        if not (len(r["ref"]) == 1 and len(r["alt"]) == 1):
            continue
        # keep even if AF missing; but MERGED_AF will show "None"
        cds, codon_anchor, _, codon_ref_cds = r["ctx"]
        key = (r["chrom"], cds["gene"], cds["strand"], codon_anchor, codon_ref_cds)
        groups[key].append(i)

    to_skip = set()
    merged = {}  # idx -> merged_text (may contain multiple lines)

    for key, idxs in groups.items():
        if len(idxs) not in (2, 3):
            continue

        chrom, gene, strand, codon_anchor, codon_ref_cds = key

                # Build multi-allelic-aware sites (unique positions) with CDS frame
        by_pos = {}
        for i in idxs:
            r = recs[i]
            frame = r["ctx"][2]
            pos = r["pos"]
            refb = r["ref"].upper()
            altb = r["alt"].upper()
            if pos not in by_pos:
                by_pos[pos] = {"pos": pos, "frame": frame, "alleles": set([refb])}
            by_pos[pos]["alleles"].add(altb)

        sites = sorted(by_pos.values(), key=lambda x: x["pos"])

        phase_res = phase_haplotypes_sites_pileup(
            bam_path=bam_path,
            chrom=chrom,
            sites=sites,
            min_mapq=min_mapq,
            min_baseq=min_baseq,
            min_informative_reads=min_informative_reads,
        )

        hap_counts = phase_res["hap_counts"]
        informative = phase_res["informative_reads"]
        dp_min = phase_res["dp_min"]
        positions = [s["pos"] for s in sites]

        # Log distribution
        if informative > 0:
            items = []
            for pat, ct in sorted(hap_counts.items(), key=lambda kv: (-kv[1], kv[0])):
                items.append(f"{pat}:{ct}:{ct/informative:.3f}")
            sys.stderr.write(
                f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
                f"[CODON_HAP] {chrom}:{min(positions)}-{max(positions)} "
                f"GENE={gene} STRAND={strand} DPmin={dp_min} INF={informative} "
                f"{' '.join(items)}\n"
            )

        else:
            sys.stderr.write(
                f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
                f"[CODON_HAP] {chrom}:{min(positions)}-{max(positions)} "
                f"GENE={gene} STRAND={strand} DPmin={dp_min} INF=0\n"
            )

        sys.stderr.flush()

                # Select haplotypes to emit (skip reference haplotype by comparing CODON_ALT to CODON_REF)
        hap_patterns = []
        if informative > 0 and phase_res.get("ok_min_reads"):
            for hap, ct in hap_counts.items():
                af_hap = ct / informative
                if af_hap < hap_af_min:
                    continue
                codon_alt_cds = codon_alt_from_bases(codon_ref_cds, strand, sites, hap)
                if codon_alt_cds == codon_ref_cds:
                    continue  # don't emit reference haplotype
                hap_patterns.append((hap, int(ct), float(af_hap), codon_alt_cds))

        if not hap_patterns:
            sys.stderr.write(
                f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] "
                f"[CODON_HAP] {chrom}:{min(positions)}-{max(positions)} "
                f"GENE={gene} STRAND={strand} no_haplotype_ge_{hap_af_min:.3g} -> keep_original_snvs\n"
            )

            sys.stderr.flush()
            continue  # keep original SNVs

        # Remove originals, replace with haplotype records
        for i in idxs:
            to_skip.add(i)

        tpl_i = min(idxs, key=lambda j: ((recs[j]["af"] if recs[j]["af"] is not None else 1e9), recs[j]["pos"]))
        tpl = recs[tpl_i]
        tpl_fields_base = tpl["fields"][:]
        cds = tpl["ctx"][0]

        # Genomic REF codon and MNV position
        if strand == "+":
            mnv_pos = codon_anchor
            ref_mnv = codon_ref_cds
            idx0 = (codon_anchor - cds["start"]) - cds["phase"]
        else:
            mnv_pos = codon_anchor - 2
            ref_mnv = revcomp(codon_ref_cds)
            idx0 = (cds["end"] - codon_anchor) - cds["phase"]
        aa_pos = (idx0 // 3) + 1

        merged_pos_str = ",".join(str(s["pos"]) for s in sites)

        # MERGED_AF: per-position AFs; if multiple ALTs exist at same POS, join AFs with '|'
        pos_to_afs = defaultdict(list)
        for i in idxs:
            pos_to_afs[recs[i]["pos"]].append(recs[i]["af"])
        merged_af_parts = []
        for s in sites:
            afs = pos_to_afs.get(s["pos"], [])
            merged_af_parts.append("|".join(str(a) for a in afs))
        merged_af_str = ",".join(merged_af_parts)

        out_lines = []
        for hap, ct, af_hap, codon_alt_cds in sorted(hap_patterns, key=lambda x: x[2], reverse=True):
            alt_mnv = codon_alt_cds if strand == "+" else revcomp(codon_alt_cds)

            tpl_fields = tpl_fields_base[:]
            tpl_fields[1] = str(mnv_pos)
            tpl_fields[3] = ref_mnv
            tpl_fields[4] = alt_mnv

            info_d = dict(tpl["info_d"])
            info_d["AF"] = f"{af_hap:.6g}"
            info_d["DP_CODON"] = str(dp_min)
            info_d["MERGED_POS"] = merged_pos_str
            info_d["MERGED_AF"] = merged_af_str
            info_d["CODON_REF"] = codon_ref_cds
            info_d["CODON_ALT"] = codon_alt_cds
            info_d["AA_REF"] = aa(codon_ref_cds)
            info_d["AA_ALT"] = aa(codon_alt_cds)
            info_d["GENE"] = gene
            info_d["HAP"] = hap
            info_d["HAP_CT"] = str(ct)
            info_d["HAP_DP"] = str(informative)

            new_cons = consequence_from_aa(info_d["AA_REF"], info_d["AA_ALT"])
            nuc_change = f"{mnv_pos}{ref_mnv}>{alt_mnv}"

            if "BCSQ" in info_d:
                info_d["BCSQ"] = update_bcsq_value(
                    info_d["BCSQ"], new_cons, aa_pos,
                    info_d["AA_REF"], info_d["AA_ALT"], nuc_change
                )

            tpl_fields[7] = format_info(info_d)
            out_lines.append("\t".join(tpl_fields))

        merged[min(idxs)] = "\n".join(out_lines)

    # Output
    header = ensure_info_headers(header)
    for h in header:
        print(h)

    for i, r in enumerate(recs):
        if i in merged:
            sys.stdout.write(merged[i] + ("\n" if not merged[i].endswith("\n") else ""))
        elif i in to_skip:
            continue
        else:
            print(r["line"])

if __name__ == "__main__":
    if len(sys.argv) < 5:
        sys.stderr.write(
            "Usage: merge_codon_mutations.py <in.vcf[.gz] or -> <ref.fasta> <ref.gff> <bam> [hap_af_min] [min_informative_reads] [min_mapq] [min_baseq]\n"
            "  - BAM is REQUIRED.\n"
            "  - hap_af_min: emit non-ref haplotypes with freq >= hap_af_min (default 0.10).\n"
            "  - min_informative_reads: minimum reads covering all positions to consider codon (default 10).\n"
            "  - min_mapq: minimum read mapping quality (default 40).\n"
            "  - min_baseq: minimum base quality at each SNV position (default 20).\n"
        )
        sys.exit(2)

    vcf_path, fasta_path, gff_path, bam_path = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
    hap_af_min = float(sys.argv[5]) if len(sys.argv) >= 6 else 0.10
    min_informative_reads = int(sys.argv[6]) if len(sys.argv) >= 7 else 10
    min_mapq = int(sys.argv[7]) if len(sys.argv) >= 8 else 40
    min_baseq = int(sys.argv[8]) if len(sys.argv) >= 9 else 20

    main(
        vcf_path, fasta_path, gff_path, bam_path,
        hap_af_min=hap_af_min,
        min_informative_reads=min_informative_reads,
        min_mapq=min_mapq,
        min_baseq=min_baseq,
    )