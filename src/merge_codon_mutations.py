#!/usr/bin/env python3
"""
merge_codon_mutations.py
- Fusionne les SNV (AF>=seuil) qui tombent dans le même codon en un seul variant
- Recalcule CODON_REF/CODON_ALT et AA_REF/AA_ALT
- La ligne fusionnée reprend comme "template" la ligne du variant le moins fréquent (AF minimal).

"""
import sys, gzip
from collections import defaultdict

# pysam est requis uniquement pour la vérification de linkage via BAM.
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
        except:
            pass
    # FORMAT/AF
    if fmt_keys and sample_fields and "AF" in fmt_keys:
        try:
            return float(sample_fields[fmt_keys.index("AF")].split(",")[0])
        except:
            pass
    # FORMAT/AD
    if fmt_keys and sample_fields and "AD" in fmt_keys:
        try:
            ad = sample_fields[fmt_keys.index("AD")].split(",")
            if len(ad) >= 2:
                r = float(ad[0]); a = float(ad[1])
                return a / (r + a) if (r + a) > 0 else None
        except:
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

def phase_linkage_snvs(
    bam_path: str,
    chrom: str,
    snvs,
    min_mapq: int = 40,
    min_baseq: int = 20,
    min_informative_reads: int = 10,
):
    """Vérifie si 2-3 SNV sont majoritairement portés ensemble sur les mêmes reads.

    snvs: liste de dicts {pos (1-based), ref, alt} en coordonnées génome (VCF).
    Linkage = AAA/(AAA+MIXED) où:
      - AAA: tous ALT sur les positions
      - MIXED: au moins un ALT mais pas tous ALT
    """
    if pysam is None:
        raise RuntimeError("pysam n'est pas disponible mais bam_path a été fourni (installez pysam).")

    pos0s = [s["pos"] - 1 for s in snvs]
    start0 = min(pos0s)
    end0 = max(pos0s)

    bam = pysam.AlignmentFile(bam_path, "rb")
    counts = defaultdict(int)
    informative = 0

    for read in bam.fetch(chrom, start0, end0 + 1):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        if read.query_qualities is None:
            continue

        states = []
        ok = True
        for s_ in snvs:
            base, qpos = _read_base_at_refpos(read, s_["pos"] - 1)
            if base is None or qpos is None:
                ok = False
                break
            if read.query_qualities[qpos] < min_baseq:
                ok = False
                break
            b = base.upper()
            if b not in "ACGT":
                ok = False
                break

            refb = s_["ref"].upper()
            altb = s_["alt"].upper()
            if b == refb:
                states.append("R")
            elif b == altb:
                states.append("A")
            else:
                states.append("O")

        if not ok:
            continue

        if "O" in states:
            counts["OTHER"] += 1
            informative += 1
            continue

        if all(x == "R" for x in states):
            counts["RRR"] += 1
        elif all(x == "A" for x in states):
            counts["AAA"] += 1
        else:
            if "A" in states:
                counts["MIXED"] += 1
            else:
                counts["OTHER"] += 1
        informative += 1

    bam.close()

    alt_support = counts["AAA"] + counts["MIXED"]
    linkage = (counts["AAA"] / alt_support) if alt_support > 0 else 0.0

    return {
        "counts": dict(counts),
        "informative_reads": informative,
        "alt_support_reads": alt_support,
        "linkage_AAA_over_alt": linkage,
        "ok_min_reads": informative >= min_informative_reads,
    }

def should_merge_minor_mnv(phase_res: dict, linkage_thr: float = 0.90) -> bool:
    return bool(phase_res.get("ok_min_reads")) and float(phase_res.get("linkage_AAA_over_alt", 0.0)) >= linkage_thr

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
        ref_b = revcomp(ref_base.upper())
        alt_b = revcomp(alt_base.upper())
    else:
        ref_b = ref_base.upper()
        alt_b = alt_base.upper()
    # on applique même si mismatch ref (VCF pas aligné sur même ref)
    cod[frame] = alt_b
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

def main(vcf_path, fasta_path, gff_path, af_thr=0.5, bam_path=None, linkage_thr=0.9, min_informative_reads=10, min_mapq=40, min_baseq=20):
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

        # Group SNVs by codon (inclut aussi les minoritaires si bam_path est fourni)
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
    merged = {}  # idx -> merged_line (print at earliest idx)

    for key, idxs in groups.items():
        # On fusionne seulement si 2 ou 3 SNV dans le même codon
        if len(idxs) not in (2, 3):
            continue

        # Décision merge:
        # - si toutes les AF >= af_thr => merge direct (historique)
        # - sinon, si bam_path fourni => merge seulement si linkage sur reads >= linkage_thr
        afs = [recs[i]["af"] for i in idxs]
        all_major = all(a is not None and a >= af_thr for a in afs)

        do_merge = False
        if all_major:
            do_merge = True
        elif bam_path:
            chrom = recs[idxs[0]]["chrom"]
            snvs = [{"pos": recs[i]["pos"], "ref": recs[i]["ref"], "alt": recs[i]["alt"]} for i in idxs]
            phase_res = phase_linkage_snvs(
                bam_path=bam_path,
                chrom=chrom,
                snvs=snvs,
                min_mapq=min_mapq,
                min_baseq=min_baseq,
                min_informative_reads=min_informative_reads,
            )

            # LOG linkage pour ce codon testé
            positions = [s["pos"] for s in snvs]
            aaa = phase_res["counts"].get("AAA", 0)
            mixed = phase_res["counts"].get("MIXED", 0)
            inform = phase_res["informative_reads"]
            alt_support = phase_res["alt_support_reads"]
            linkage = phase_res["linkage_AAA_over_alt"]

            sys.stderr.write(
                f"[PHASE] {chrom}:{min(positions)}-{max(positions)} "
                f"GENE={key[1]} STRAND={key[2]} "
                f"AAA={aaa} MIXED={mixed} ALT={alt_support} INF={inform} "
                f"LINK={linkage:.3f}\n"
            )
            sys.stderr.flush()

            do_merge = should_merge_minor_mnv(phase_res, linkage_thr=linkage_thr)

        if not do_merge:
            continue
        chrom, gene, strand, codon_anchor, codon_ref_cds = key

        # construire codon_alt
        codon_alt_cds = codon_ref_cds
        pos_list = []
        af_list = []
        for i in idxs:
            r = recs[i]
            cds, _, frame, _ = r["ctx"]
            codon_alt_cds = apply_snv(codon_alt_cds, frame, r["ref"], r["alt"], strand)
            pos_list.append(r["pos"])
            af_list.append(r["af"])

        # template = variant le moins fréquent (AF min)
        tpl_i = min(idxs, key=lambda j: (recs[j]["af"], recs[j]["pos"]))
        tpl = recs[tpl_i]
        tpl_fields = tpl["fields"][:]

        cds = tpl["ctx"][0]

        # construire REF/ALT génomiques sur 3bp (toujours codon complet)
        if strand == "+":
            mnv_pos = codon_anchor
            ref_mnv = codon_ref_cds
            alt_mnv = codon_alt_cds
            codon_start = codon_anchor
            idx0 = (codon_start - cds["start"]) - cds["phase"]
        else:
            mnv_pos = codon_anchor - 2
            ref_mnv = revcomp(codon_ref_cds)
            alt_mnv = revcomp(codon_alt_cds)
            idx0 = (cds["end"] - codon_anchor) - cds["phase"]
        aa_pos = (idx0 // 3) + 1
        tpl_fields[1] = str(mnv_pos)
        tpl_fields[3] = ref_mnv
        tpl_fields[4] = alt_mnv

        info_d = dict(tpl["info_d"])  # garde le moins fréquent
        info_d["AF"] = str(min(af_list))  # conservateur
        info_d["MERGED_POS"] = ",".join(map(str, sorted(pos_list)))
        # garder MERGED_AF dans le même ordre que MERGED_POS
        pos_to_af = {recs[i]["pos"]: recs[i]["af"] for i in idxs}
        info_d["MERGED_AF"] = ",".join(str(pos_to_af[p]) for p in sorted(pos_list))
        info_d["CODON_REF"] = codon_ref_cds
        info_d["CODON_ALT"] = codon_alt_cds
        info_d["AA_REF"] = aa(codon_ref_cds)
        info_d["AA_ALT"] = aa(codon_alt_cds)
        info_d["GENE"] = gene

        # Recalcul consequence + réécriture BCSQ si présent
        new_cons = consequence_from_aa(info_d["AA_REF"], info_d["AA_ALT"])

        # nuc change codon-level en coordonnées génome
        # ex: "898AAT>GGC" (utilise mnv_pos/ref_mnv/alt_mnv déjà calculés)
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
            print(merged[i])
        elif i in to_skip:
            continue
        else:
            print(r["line"])

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

if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.stderr.write(
            "Usage: merge_codon_mutations.py <in.vcf[.gz] or -> <ref.fasta> <ref.gff> [af_threshold] [bam_path] [linkage_threshold]\n"
            "  - af_threshold: seuil 'majoritaire' (défaut 0.5)\n"
            "  - bam_path: si fourni, permet de fusionner des MNV minoritaires après test de linkage sur reads\n"
            "  - linkage_threshold: seuil AAA/(AAA+MIXED) (défaut 0.9)\n"
        )
        sys.exit(2)

    vcf_path, fasta_path, gff_path = sys.argv[1], sys.argv[2], sys.argv[3]
    af_thr = float(sys.argv[4]) if len(sys.argv) >= 5 and sys.argv[4] not in ("-", "") else 0.5
    bam_path = sys.argv[5] if len(sys.argv) >= 6 and sys.argv[5] not in ("-", "") else None
    linkage_thr = float(sys.argv[6]) if len(sys.argv) >= 7 else 0.9

    main(
        vcf_path, fasta_path, gff_path,
        af_thr=af_thr,
        bam_path=bam_path,
        linkage_thr=linkage_thr,
    )
