# hgene tests

This directory contains all functional tests for the **hgene** pipeline.

Each test validates a specific script or logical component of the
workflow.\
All datasets are synthetic and deterministic.

------------------------------------------------------------------------

# 1. Mapping & Core Pipeline

------------------------------------------------------------------------

## test_01_map2HHV.sh

**Tests:** `hgene_map.sh`\
Basic validation of read mapping against the HHV1 reference.

------------------------------------------------------------------------

## test_02_map2HHV.sh

**Tests:** `hgene_map.sh`\
Checks mapping reproducibility and output consistency.

------------------------------------------------------------------------

## test_03_expected_vcf.sh

**Tests:** `hgene_variant_call.sh` (VCF production)\
Checks that the produced VCF matches the expected structure/content for
the test dataset.

------------------------------------------------------------------------

# 2. Custom Filtering Logic

------------------------------------------------------------------------

These tests validate `hgene_filter_custom.py`.

## test_04_hpoly_deletion_removed.sh

Deletion in homopolymer context removed when AF below threshold.

------------------------------------------------------------------------

## test_05_noHpoly_deletion_removed.sh

Deletion removed when not in homopolymer and AF below threshold.

------------------------------------------------------------------------

## test_06_noHpoly_deletion_kept.sh

Deletion retained when AF above threshold no homopolymer (HRUN=3).

------------------------------------------------------------------------

## test_07_noHpoly_deletion_kept.sh

Deletion retained when AF above threshold no homopolymer (no HRUN).


------------------------------------------------------------------------

## test_08_hpoly_deletion_kept.sh

Deletion retained in homopolymer when AF satisfies condition.

------------------------------------------------------------------------

## test_09_minVar_removed.sh

Variant removed when AF \< 0.5 and depth below required minimum.

------------------------------------------------------------------------

## test_10_minVar_kept.sh

Variant retained when AF \< 0.5 but depth requirement is satisfied.

------------------------------------------------------------------------

# 3. Codon Haplotype Merge Logic

These tests validate `hgene_codon_haplotype_merge.py`.

Merging is based on: **Read-level phasing**.

------------------------------------------------------------------------

## test_11_merge_two_mut_same_codon.sh

Scenario: Two SNVs within the same codon.

Expected: - Single codon-level variant emitted
- Original SNVs removed
- HAP, HAP_CT, HAP_DP present

------------------------------------------------------------------------

## test_12_merge_three_mut_same_codon.sh

Scenario: Three SNVs in the same codon.

Expected: - Fully merged codon-level variant

------------------------------------------------------------------------

## test_13_no_merge_diff_gene.sh

Scenario: Two SNVs located in different genes (but same codon position).

Expected: - No merge

------------------------------------------------------------------------

## test_14_no_merge_adjacent_codons.sh

Scenario: Two SNVs in adjacent codons.

Expected: - No merge

------------------------------------------------------------------------

# 4. Advanced Population Scenarios

These tests validate complex haplotype reconstruction.

------------------------------------------------------------------------

## test_15_mixed_linked_unlinked_populations.sh

Scenario: Mixed linked and unlinked populations.

### Expected codon distribution

| Type | Count | AF  | Codon |
|------|-------|-----|-------|
| REF  | 40    | 0.3 | TCG   |
| ALT  | 20    | 0.2 | GAG   |
| ALT  | 30    | 0.3 | GCG   |
| ALT  | 10    | 0.1 | TAG   |


Output behavior: - Only non-reference codons emitted
- No collapsing of populations
- Correct HAP_CT and HAP_DP values

------------------------------------------------------------------------

## test_16_multiallelic_same_pos.sh

Scenario: Multi-allelic variants at the same genomic position.

Expected codon distribution:

### Expected codon distribution

| Type | Count | AF  | Codon |
|------|-------|-----|-------|
| REF  | 30    | 0.3 | TCG   |
| ALT  | 20    | 0.2 | GAG   |
| ALT  | 30    | 0.3 | GCG   |
| ALT  | 10    | 0.1 | TAG   |
| ALT  | 10    | 0.1 | ACG   |


Output behavior: - All haplotypes ≥ threshold emitted
- Reference codon excluded
- Multi-allelic support verified
- No collapsing of populations

------------------------------------------------------------------------

