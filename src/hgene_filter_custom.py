#!/usr/bin/env python3
import sys

THR_DEFAULT = 0.20       # deletion si HRUN < 4 (ou absent)
THR_HRUN_4PLUS = 0.40    # deletion si HRUN >= 4
MIN_DP_MINOR = 100       # si AF < 0.5, DP minimal

def get_info(info, key):
    for f in info.split(";"):
        if f.startswith(key + "="):
            return f.split("=", 1)[1]
    return None

vcf_path = sys.argv[1]
printed_filters = False

with open(vcf_path, "r") as vcf:
    for l in vcf:
        if l.startswith("##"):
            print(l.rstrip())
            continue
        # Column header line: inject filter definitions once, just before #CHROM
        if l.startswith("#CHROM"):
          if not printed_filters:
                    print(
                        f'##FILTER=<ID=min_dp_minor,Description="Minor variants (AF < 0.5) require DP >= {MIN_DP_MINOR}">'
                    )
                    print(
                        f'##FILTER=<ID=del_hrun_ge4,Description="Deletions require AF >= {THR_HRUN_4PLUS} when HRUN >= 4">'
                    )
                    print(
                        f'##FILTER=<ID=del_hrun_lt4,Description="Deletions require AF >= {THR_DEFAULT} when HRUN < 4 or HRUN missing">'
                    )
                    printed_filters = True
          print(l.rstrip())
          continue
        if l.startswith("#"):
            print(l.rstrip())
            continue


        ligne = l.rstrip("\n").split("\t")
        ref = ligne[3]
        alt = ligne[4]
        info = ligne[7]

        af_s = get_info(info, "AF")
        dp_s = get_info(info, "DP")
        hrun_s = get_info(info, "HRUN")

        # Si AF absent, on ne peut pas appliquer les filtres -> skip
        if af_s is None:
            continue

        af = float(af_s)
        dp = int(dp_s) if dp_s is not None else 0
        hrun = int(hrun_s) if hrun_s is not None else 0

        # variants minoritaires
        if af < 0.5 and dp < MIN_DP_MINOR:
            continue

        # Si pas deletion -> print
        if len(ref) <= len(alt):
            print(l.rstrip())
            continue

        # Deletion: seuil dépend de HRUN
        thr = THR_HRUN_4PLUS if hrun >= 4 else THR_DEFAULT
        if af >= thr:
            print(l.rstrip())
