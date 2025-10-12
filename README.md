# fastq-filter — bioinformatics utilities (DNA/RNA + FASTQ filtering)

Small educational toolkit for working with nucleotide sequences and FASTQ files.  

- DNA/RNA tools: reverse, complement, reverse_complement, transcribe, is_nucleic_acid
- FASTQ filtering by GC%, read length, and mean quality (Phred+33)
- File-based streaming mode — processes large FASTQ files without loading them into memory

> The entry file exposes exactly two main functions:  
> run_dna_rna_tools(...) and filter_fastq(...).

---

## Repository structure
.
├── fastq_filter.py             # Entry point (CLI + 2 facades)
└── modules/
    ├── dna_rna_tools_module.py # DNA/RNA helper functions (pure)
    ├── fastq_utils.py          # FASTQ metrics and filtering logic
    ├── in_out_utils.py         # File I/O helpers for FASTQ reading/writing
    ├── example_data.py
    └── example_fastq.fastq

---

## Installation / Development setup
git clone https://github.com/alinatgrv/fastq-filter.git
cd fastq-filter

---

## DNA/RNA tools

from fastq_filter import run_dna_rna_tools

run_dna_rna_tools("ATGC", "reverse_complement")   # "GCAT"
run_dna_rna_tools("AUGC", "transcribe")           # "ATGC"
run_dna_rna_tools("ATGC", "is_nucleic_acid")      # True
run_dna_rna_tools("ATGC", "GGG", "reverse")       # ["CGTA", "GGG"]

Supported procedures:
is_nucleic_acid | transcribe | reverse | complement | reverse_complement

---

## FASTQ Filtering

Filters FASTQ (or FASTQ.GZ) directly from disk using streaming iteration.

### Example
python3 fastq_filter.py \
  --input_fastq modules/example_fastq.fastq \
  --output_fastq kept.fastq \
  --gc_bounds 40 60 \
  --length_bounds 50 300 \
  --quality_threshold 20

Output (example):
Kept 89 / 89 reads. Saved to: filtered/kept.fastq

### Command-line arguments
Argument | Type / Default | Description
--------- | --------------- | ------------
--input_fastq | required | Path to input FASTQ or FASTQ.GZ file
--output_fastq | required | Output file name (created under `filtered/`)
--gc_bounds | 0 100 | GC% bounds (one value = upper bound; two = lower upper)
--length_bounds | 0 4294967296 | Read length bounds (inclusive)
--quality_threshold | 0.0 | Minimum mean Phred+33 quality
--force | False | Allow overwriting existing file (otherwise _1, _2, … suffix)

### Output directory
Filtered reads are saved to the folder filtered/.  
If the output file already exists, a new file with a numeric suffix is created.

Example:
filtered/
├── kept.fastq
├── kept_1.fastq
└── kept_2.fastq

### Implementation details
- Streamed reading and writing handled by modules/in_out_utils.py
- Filtering logic (GC%, length, quality) implemented in modules/fastq_utils.py
- Entry CLI logic located in fastq_filter.py
- Works both with .fastq and .fastq.gz

---

## API Summary

Function | Description
--------- | ------------
run_dna_rna_tools(*args: str) | Applies the requested DNA/RNA transformation. Last argument = method name.
filter_fastq(input_fastq_path, output_fastq_name, ...) | Reads a FASTQ file, filters reads on-the-fly, and writes results to filtered/.

---

## References
- Illumina Quality Score Encoding (Phred+33): https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/QualityScoreEncoding_swBS.htm  
- Python ord() documentation: https://docs-python.ru/tutorial/vstroennye-funktsii-interpretatora-python/funktsija-ord/  
- ASCII table: https://www.asciitable.com/

---

## Author & Contacts

Alina Tagirova  
📧 hortizaira@gmail.com  
💻 GitHub: https://github.com/alinatgrv

---

This project was developed as part of the Institute of Bioinformatics course (2025).  
It demonstrates practical skills in Python, modular design, and working with biological file formats.
