# hgene

**Version:** 1.0.1-rc1

hgene performs small-variant analysis (SNPs, MNVs and indels) in the HSV resistance genes UL23 and UL30, applies homopolymer-aware filtering, reconstructs codon-level amino-acid consequences using read-level linkage, and reports variant co-occurrence evidence.

---

## Usage

`hgene -v <virus> [-c <cpu>] <prefix|prefix.fastq>`

**Arguments**

- `-v <virus>` — Virus reference key (e.g. HHV1, HHV2)
- `-c <cpu>` — Number of threads (default: nproc)
- `<input>` — Prefix or an uncompressed `.fastq` file

**Notes**

- `.fq` files are not supported
- Compressed files (`.gz`) are not supported

---

## Key parameters and filtering thresholds

### Read preprocessing
- Adapter trimming performed with `porechop --discard_middle`  
_- Reads containing internal adapters are fully discarded_

### Mapping
- `minimap2 -x map-ont`  

### BAM filtering
- Minimum mapping quality: **MAPQ ≥ 40**
- Minimum aligned reference span: **≥ 1000 bp**

### Variant calling (LoFreq)
- Default LoFreq filters disabled
- Indel calling enabled
- Initial filtering:
  - Minimum coverage: **20**
  - Minimum allele frequency (AF): **0.10**

### Custom filtering (homopolymers and minor variants)
- Minor variants (`AF < 0.5`) require **DP ≥ 100**
- Deletions:
  - **AF ≥ 0.40** when `HRUN ≥ 4`
  - **AF ≥ 0.20** when `HRUN < 4` or missing

### Codon reconstruction and linkage
- Codon-level variant reconstruction uses read-level evidence
- Minimum haplotype allele frequency: **0.10**
- Minimum informative reads: **10**
- Minimum MAPQ for linkage: **40**
- Minimum base quality: **20**

### Output indexing
- Final VCF compressed and indexed using **CSI**

---

## Output structure

A run produces a single directory named:

```
sample_output/
```

### Directory tree

```
sample_output/
├── BAM_sample
│   ├── sample.bam
│   └── sample.bam.bai
├── REPORT_sample
│   ├── sample_bam.html
│   ├── sample.html
│   └── sample.vcf.gz.<virus>.docx
├── VCF_sample
│   ├── sample.bcf.vcf.gz
│   ├── sample.filtered.vcf.gz
│   ├── sample.filtered.vcf.gz.csi
│   ├── sample.lofreq_filtered.vcf
│   └── sample.lofreq_raw.vcf
├── sample.vcf.gz
└── sample.vcf.gz.csi
```

### Contents

**Root**
- `sample.vcf.gz` : final result  
- `sample.vcf.gz.csi` : index  

**REPORT_sample/**     
_In french_ 
- `sample_bam.html` : mapping report  
- `sample.html` : variant report  
- `sample.vcf.gz.<virus>.docx` : clinical report  

**BAM_sample/**
- `sample.bam` : BAM file  
- `sample.bam.bai` : BAM index  

**VCF_sample/**  
 _All intermediate VCF files_

- `sample.lofreq_raw.vcf`
- `sample.lofreq_filtered.vcf`
- `sample.filtered.vcf.gz`
- `sample.filtered.vcf.gz.csi`
- `sample.bcf.vcf.gz`

---

## Custom report template

Users can provide a custom `template.docx` by placing it in `$HOME/template/hsv-1/` or `$HOME/template/hsv-2/`. If detected, it will automatically be used for report generation. If no template is found, the default built-in template is used.

---

## Acknowledgments

Julien Prados (@pradosj, Bioinformatics Support Platform) — for his contribution.

---

## Disclaimer

```
IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
```
