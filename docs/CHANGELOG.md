# Changelog

All notable changes to **hgene** are documented in this file.

## [1.1.0]

### Added
- Support for CMV reporting, including HHV5 genes and resistance interpretation.
- Parallelized Porechop preprocessing using the configured CPU/thread count.

### Changed
- Merged the previous HSV and CMV branches into a single codebase.
- Updated pipeline structure to use external database folders instead of hardcoded internal references.
- Improved traceability of analyses through database and analysis metadata.

---

## [1.0.1]

### Added
- Integration of virotyper for variant interpretation and reporting.
- Codon-level multi-nucleotide variant reconstruction using read-level linkage.
- Metadata reporting in VCF header.

### Changed
- Improved minor variant filtering requiring DP ≥ 100 when AF < 0.5
- Refined deletion allele-frequency thresholds (HRUN-dependent)
- Improved README documentation and parameter transparency

---

## [1.0.0] — Initial release
