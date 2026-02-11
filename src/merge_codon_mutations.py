#!/usr/bin/env python3
"""
merge_codon_mutations.py

- Mode majoritaire seulement:
    - Fusionne les SNV (AF>=seuil) qui tombent dans le même codon en un seul variant
    - Recalcule CODON_REF/CODON_ALT et AA_REF/AA_ALT
    - La ligne fusionnée reprend comme "template" la ligne du variant le moins fréquent (AF minimal).

- Mode haplotypes codon (initialement pour minoritaire) :
    - Pour chaque codon contenant 2 ou 3 SNV, reconstruit les haplotypes read-level (patterns R/A)
    - Reporte TOUTE combinaison (≠ ref) dont la fréquence >= hap_af_min (défaut 0.10)
    - AF est recalculée comme: count(haplotype)/informative_reads
    - DP_CODON est min(depth_position) sur les positions du codon (après filtres MAPQ/baseQ)
    - Ajoute HAP/HAP_CT/HAP_DP dans INFO.
    - HAP_DP est le nombre de reads couvrant toutes les positions du codon (profondeur réellement utilisée pour calculer l’AF haplotypique), tandis que DP_CODON est la profondeur minimale par position dans le codon, donc théoriquement ≥ HAP_DP
    - DP devrait être remplacé par HAP_DP ?
    - Log en stderr: distribution des haplotypes et fréquences.

Usage:
  merge_codon_mutations.py <in.vcf[.gz] or -> <ref.fasta> <ref.gff> [af_threshold] [bam_path] [linkage_threshold] [hap_af_min]

Notes:
  - linkage_threshold utilisé seulement si hap_af_min est None)
  - hap_af_min (float) active le mode haplotypes, si bam_path est fourni (défaut 0.10)
"""
import sys, gzip
from collections import defaultdict

# pysam est requis pour la vérification via BAM.
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
    """Retourne (base, qpos) au ref_pos0 (0-based), ou (None, None) si non couvert."""
    qpos_at = None
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if rpos == ref_pos0:
            if qpos is None:
                return None, None
            qpos_at = qpos
            break
    if qpos_at is None:
        return None, None
    return read.query_sequence[qpos_at], qpos_at

def phase_haplotypes_snvs(
    bam_path: str,
    chrom: str,
    snvs,
    min_mapq: int = 40,
    min_baseq: int = 20,
    min_informative_reads: int = 10,
):
    """
    Reconstruit les haplotypes read-level sur 2-3 SNV.
    Retourne:
      - hap_counts: patterns (ex "AR", "AA", "RRR", ...) + "OTHER"
      - informative_reads: nb reads couvrant toutes les positions avec base A/C/G/T et baseQ>=min_baseq
      - dp_pos: profondeur "valide" par position (A/C/G/T + baseQ>=min_baseq) parmi reads MAPQ>=min_mapq
      - dp_min: min(dp_pos)
    """
    if pysam is None:
        raise RuntimeError("pysam is not available but bam_path has been provided (install pysam).")

    pos0s = [s["pos"] - 1 for s in snvs]
    start0 = min(pos0s)
    end0 = max(pos0s)

    bam = pysam.AlignmentFile(bam_path, "rb")
    hap_counts = defaultdict(int)
    dp_pos = [0] * len(snvs)
    informative = 0

    for read in bam.fetch(chrom, start0, end0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.query_qualities is None:
            continue

        states = []
        ok_all_positions = True

        for j, s_ in enumerate(snvs):
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

            # DP par position: read couvre cette position avec base "propre"
            dp_pos[j] += 1

            refb = s_["ref"].upper()
            altb = s_["alt"].upper()
            if b == refb:
                states.append("R")
            elif b == altb:
                states.append("A")
            else:
                states.append("O")

        if not ok_all_positions:
            continue

        informative += 1

        if "O" in states:
            hap_counts["OTHER"] += 1
        else:
            hap_counts["".join(states)] += 1

    bam.close()

    dp_min = min(dp_pos) if dp_pos else 0

    return {
        "hap_counts": dict(hap_counts),
        "informative_reads": informative,
        "dp_pos": dp_pos,
        "dp_min": dp_min,
        "ok_min_reads": informative >= min_informative_reads,
    }

def should_merge_minor_mnv_linkage(phase_res: dict, linkage_thr: float = 0.90) -> bool:
    """Compat: ancien critère AAA/(AAA+MIXED) à partir de hap_counts."""
    if not phase_res.get("ok_min_reads"):
        return False
    hap_counts = phase_res.get("hap_counts", {})
    # Construire AAA/MIXED pour 2-3 positions
    # AAA: pattern tout A
    # MIXED: patterns avec au moins un A mais pas tous A
    patterns = [p for p in hap_counts.keys() if p not in ("OTHER",)]
    if not patterns:
        return False
    # déduire longueur
    L = None
    for p in patterns:
        if p != "OTHER":
            L = len(p)
            break
    if L is None:
        return False
    aaa_pat = "A" * L
    aaa = int(hap_counts.get(aaa_pat, 0))
    mixed = 0
    for p, ct in hap_counts.items():
        if p in ("OTHER", aaa_pat):
            continue
        if "A" in p and p != ("R"*L):
            mixed += int(ct)
    alt_support = aaa + mixed
    linkage = (aaa / alt_support) if alt_support > 0 else 0.0
    return linkage >= linkage_thr

def codon_ctx(chrom, pos1, fasta, cds_by_chr):
    """
    Retourne (cds, codon_anchor, frame, codon_ref_CDS)
    - pos1: position 1-based génomique
    - codon_anchor: + => start du codon (1-based, gauche)
                    - => position la plus à droite du triplet sur le génome (1-based)
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
                codon_anchor = pos1 + frame  # rightmost base of triplet on genome
                triplet = seq[codon_anchor-3:codon_anchor]
                if len(triplet) != 3:
                    return None
                codon = revcomp(triplet)
                return (cds, codon_anchor, frame, codon)
    return None

def apply_snv(codon_ref_cds, frame, ref_base, alt_base, strand):
    cod = list(codon_ref_cds)
    if strand == "-":
        alt_b = revcomp(alt_base.upper())
    else:
        alt_b = alt_base.upper()
    cod[frame] = alt_b
    return "".join(cod)

def codon_alt_from_pattern(codon_ref_cds: str, strand: str, snvs_with_frame, pattern: str) -> str:
    """Construit CODON_ALT_CDS en appliquant uniquement les SNV où pattern a 'A'."""
    cod = list(codon_ref_cds)
    for s, st in zip(snvs_with_frame, pattern):
        if st != "A":
            continue
        alt_base = s["alt"].upper()
        if strand == "-":
            alt_base = revcomp(alt_base)
        cod[s["frame"]] = alt_base
    return "".join(cod)

def ensure_info_headers(header_lines):
    def has(id_):
        return any(h.startswith("##INFO=") and f"ID={id_}," in h for h in header_lines)

    extras = []
    defs = [
        ("MERGED_POS", '##INFO=<ID=MERGED_POS,Number=1,Type=String,Description="Original SNV positions merged into one codon variant">'),
        ("MERGED_AF",  '##INFO=<ID=MERGED_AF,Number=1,Type=String,Description="AFs of original SNVs (same order as MERGED_POS)">'),
        ("CODON_REF",  '##INFO=<ID=CODON_REF,Number=1,Type=String,Description="Reference codon (CDS orientation)">'),
        ("CODON_ALT",  '##INFO=<ID=CODON_ALT,Number=1,Type=String,Description="Mutated codon (CDS orientation)">'),
        ("AA_REF",     '##INFO=<ID=AA_REF,Number=1,Type=String,Description="Reference amino acid">'),
        ("AA_ALT",     '##INFO=<ID=AA_ALT,Number=1,Type=String,Description="Mutated amino acid">'),
        ("GENE",       '##INFO=<ID=GENE,Number=1,Type=String,Description="Gene/Parent/ID extracted from GFF">'),
        ("DP_CODON",   '##INFO=<ID=DP_CODON,Number=1,Type=Integer,Description="Minimum per-site depth across the merged codon positions (after MAPQ/baseQ filters)">'),
        ("HAP",        '##INFO=<ID=HAP,Number=1,Type=String,Description="Haplotype pattern across merged SNVs (R/A per site, same order as MERGED_POS)">'),
        ("HAP_CT",     '##INFO=<ID=HAP_CT,Number=1,Type=Integer,Description="Read count supporting this haplotype (pattern)">'),
        ("HAP_DP",     '##INFO=<ID=HAP_DP,Number=1,Type=Integer,Description="Informative read depth used to estimate haplotype AF">'),
    ]
    for id_, line in defs:
        if not has(id_):
            extras.append(line)

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

def update_bcsq_value(bcsq: str, new_consequence: str, strand: str, aa_pos: int,
                      aa_ref: str, aa_alt: str, nuc_change: str) -> str:
    """
    Met à jour le 1er enregistrement BCSQ.
    Format supposé:
      consequence|gene|...|protein_coding|+|300N|900T>C
    On remplace:
      [0] consequence
      [5] AA field -> "{aa_pos}{aa_ref}>{aa_pos}{aa_alt}"
      [6] nuc field -> nuc_change (ex "898AAT>GGC")
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

def main(
    vcf_path,
    fasta_path,
    gff_path,
    af_thr=0.5,
    bam_path=None,
    linkage_thr=0.9,
    hap_af_min=0.10,
    min_informative_reads=10,
    min_mapq=40,
    min_baseq=0   #permisif pour pouvoir run les test normalement la qualité est dejà filtré par le gridIon
):
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
            alt = fields[4].split(",")[0]  # simple: 1er ALT
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
                "fmt_keys": fmt_keys,
                "sample_fields": sample_fields,
                "af": af,
                "ctx": ctx,  # (cds, codon_anchor, frame, codon_ref_cds)
            })

    # Group SNVs by codon
    groups = defaultdict(list)  # key -> list of idx
    for i, r in enumerate(recs):
        if r["ctx"] is None:
            continue
        if not (len(r["ref"]) == 1 and len(r["alt"]) == 1):
            continue
        if r["af"] is None:
            continue
        cds, codon_anchor, _, codon_ref_cds = r["ctx"]
        key = (r["chrom"], cds["gene"], cds["strand"], codon_anchor, codon_ref_cds)
        groups[key].append(i)

    to_skip = set()
    merged = {}  # idx -> merged_text (can contain multiple lines)

    for key, idxs in groups.items():
        if len(idxs) not in (2, 3):
            continue

        chrom, gene, strand, codon_anchor, codon_ref_cds = key

        # MODE haplotypes (si BAM fourni)
        if bam_path:
            # Build snvs with frame in same order as positions
            snvs = []
            for i in idxs:
                r = recs[i]
                frame = r["ctx"][2]
                snvs.append({"pos": r["pos"], "ref": r["ref"], "alt": r["alt"], "frame": frame})

            # deterministic order: sort by genomic pos (and reorder frames accordingly)
            snvs = sorted(snvs, key=lambda x: x["pos"])

            phase_res = phase_haplotypes_snvs(
                bam_path=bam_path,
                chrom=chrom,
                snvs=snvs,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                min_informative_reads=min_informative_reads,
            )

            hap_counts = phase_res["hap_counts"]
            informative = phase_res["informative_reads"]
            dp_min = phase_res["dp_min"]
            ref_pattern = "R" * len(snvs)

            # Log distribution
            if informative > 0:
                # compact summary
                items = []
                for pat, ct in sorted(hap_counts.items(), key=lambda kv: (-kv[1], kv[0])):
                    if pat == "OTHER":
                        continue
                    items.append(f"{pat}:{ct}:{ct/informative:.3f}")
                sys.stderr.write(
                    f"[CODON_HAP] {chrom}:{min(s['pos'] for s in snvs)}-{max(s['pos'] for s in snvs)} "
                    f"GENE={gene} STRAND={strand} DPmin={dp_min} INF={informative} "
                    f"{' '.join(items)}\n"
                )
                sys.stderr.flush()

            # Build hap patterns to emit
            hap_patterns = []
            if informative > 0 and phase_res.get("ok_min_reads"):
                for pat, ct in hap_counts.items():
                    if pat in ("OTHER", ref_pattern):
                        sys.stderr.write(
                            f"[HAPS] {chrom}:{min(positions)}-{max(positions)} "
                            f"GENE={key[1]} STRAND={key[2]} "
                            f"no_haplotype_ge_10pct -> keep_original_snvs\n"
                        )
                        sys.stderr.flush()
                        continue
                    af_hap = ct / informative
                    if af_hap >= hap_af_min:
                        hap_patterns.append((pat, int(ct), float(af_hap)))

            if not hap_patterns:
                continue  # keep original SNVs

            # Remove original SNVs, replace with haplotype records
            for i in idxs:
                to_skip.add(i)

            # template = least frequent original SNV
            tpl_i = min(idxs, key=lambda j: (recs[j]["af"], recs[j]["pos"]))
            tpl = recs[tpl_i]
            tpl_fields_base = tpl["fields"][:]
            cds = tpl["ctx"][0]

            # Genome-level REF/ALT of codon full triplet
            if strand == "+":
                mnv_pos = codon_anchor
                ref_mnv = codon_ref_cds
                idx0 = (codon_anchor - cds["start"]) - cds["phase"]
            else:
                mnv_pos = codon_anchor - 2
                ref_mnv = revcomp(codon_ref_cds)
                idx0 = (cds["end"] - codon_anchor) - cds["phase"]
            aa_pos = (idx0 // 3) + 1

            merged_pos_str = ",".join(str(s["pos"]) for s in snvs)
            # original SNV AFs (sorted order)
            pos_to_af = {recs[i]["pos"]: recs[i]["af"] for i in idxs}
            merged_af_str = ",".join(str(pos_to_af.get(s["pos"])) for s in snvs)

            out_lines = []
            for pat, ct, af_hap in sorted(hap_patterns, key=lambda x: x[2], reverse=True):
                codon_alt_cds = codon_alt_from_pattern(codon_ref_cds, strand, snvs, pat)
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
                info_d["HAP"] = pat
                info_d["HAP_CT"] = str(ct)
                info_d["HAP_DP"] = str(informative)

                new_cons = consequence_from_aa(info_d["AA_REF"], info_d["AA_ALT"])
                nuc_change = f"{mnv_pos}{ref_mnv}>{alt_mnv}"

                if "BCSQ" in info_d:
                    info_d["BCSQ"] = update_bcsq_value(
                        info_d["BCSQ"], new_cons, strand, aa_pos,
                        info_d["AA_REF"], info_d["AA_ALT"], nuc_change
                    )

                tpl_fields[7] = format_info(info_d)
                out_lines.append("\t".join(tpl_fields))

            merged[min(idxs)] = "\n".join(out_lines)
            continue

        # MODE classique sans BAM (ou bam_path absent)
        afs = [recs[i]["af"] for i in idxs]
        all_major = all(a is not None and a >= af_thr for a in afs)
        if not all_major:
            continue

        # construire codon_alt
        codon_alt_cds = codon_ref_cds
        pos_list = []
        af_list = []
        for i in idxs:
            r = recs[i]
            frame = r["ctx"][2]
            codon_alt_cds = apply_snv(codon_alt_cds, frame, r["ref"], r["alt"], strand)
            pos_list.append(r["pos"])
            af_list.append(r["af"])

        tpl_i = min(idxs, key=lambda j: (recs[j]["af"], recs[j]["pos"]))
        tpl = recs[tpl_i]
        tpl_fields = tpl["fields"][:]
        cds = tpl["ctx"][0]

        if strand == "+":
            mnv_pos = codon_anchor
            ref_mnv = codon_ref_cds
            alt_mnv = codon_alt_cds
            idx0 = (codon_anchor - cds["start"]) - cds["phase"]
        else:
            mnv_pos = codon_anchor - 2
            ref_mnv = revcomp(codon_ref_cds)
            alt_mnv = revcomp(codon_alt_cds)
            idx0 = (cds["end"] - codon_anchor) - cds["phase"]
        aa_pos = (idx0 // 3) + 1

        tpl_fields[1] = str(mnv_pos)
        tpl_fields[3] = ref_mnv
        tpl_fields[4] = alt_mnv

        info_d = dict(tpl["info_d"])
        info_d["AF"] = str(min(af_list))
        info_d["MERGED_POS"] = ",".join(map(str, sorted(pos_list)))
        pos_to_af = {recs[i]["pos"]: recs[i]["af"] for i in idxs}
        info_d["MERGED_AF"] = ",".join(str(pos_to_af[p]) for p in sorted(pos_list))
        info_d["CODON_REF"] = codon_ref_cds
        info_d["CODON_ALT"] = codon_alt_cds
        info_d["AA_REF"] = aa(codon_ref_cds)
        info_d["AA_ALT"] = aa(codon_alt_cds)
        info_d["GENE"] = gene

        new_cons = consequence_from_aa(info_d["AA_REF"], info_d["AA_ALT"])
        nuc_change = f"{mnv_pos}{ref_mnv}>{alt_mnv}"

        if "BCSQ" in info_d:
            info_d["BCSQ"] = update_bcsq_value(
                info_d["BCSQ"], new_cons, strand, aa_pos,
                info_d["AA_REF"], info_d["AA_ALT"], nuc_change
            )

        tpl_fields[7] = format_info(info_d)
        merged_line = "\t".join(tpl_fields)

        for i in idxs:
            to_skip.add(i)
        merged[min(idxs)] = merged_line

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
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage: merge_codon_mutations.py <in.vcf[.gz] or -> <ref.fasta> <ref.gff> "
            "[af_threshold] [bam_path] [linkage_threshold] [hap_af_min]\n"
            "  - af_threshold: seuil 'majoritaire' (défaut 0.5) [mode sans BAM]\n"
            "  - bam_path: si fourni, active le mode haplotypes (AF recalculée par haplotype)\n"
            "  - linkage_threshold: conservé pour compatibilité (non utilisé en mode haplotypes)\n"
            "  - hap_af_min: seuil haplotype (défaut 0.10) en mode haplotypes\n"
        )
        sys.exit(2)

    vcf_path, fasta_path, gff_path = sys.argv[1], sys.argv[2], sys.argv[3]
    af_thr = float(sys.argv[4]) if len(sys.argv) >= 5 and sys.argv[4] not in ("-", "") else 0.5
    bam_path = sys.argv[5] if len(sys.argv) >= 6 and sys.argv[5] not in ("-", "") else None
    linkage_thr = float(sys.argv[6]) if len(sys.argv) >= 7 else 0.9
    hap_af_min = float(sys.argv[7]) if len(sys.argv) >= 8 else 0.10

    main(
        vcf_path, fasta_path, gff_path,
        af_thr=af_thr,
        bam_path=bam_path,
        linkage_thr=linkage_thr,
        hap_af_min=hap_af_min,
    )