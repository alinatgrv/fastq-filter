# fastq-filter — bioinformatics utilities (DNA/RNA + FASTQ filtering)

Small educational toolkit:
- DNA/RNA helpers (`reverse`, `complement`, `reverse_complement`, `transcribe`, validation).
- FASTQ read filtering by GC%, length, and mean quality (Phred+33).

> Entry file exposes exactly **two** functions required by the assignment:  
> `run_dna_rna_tools(...)` and `filter_fastq(...)`.  
> Pure functions only — no I/O inside modules.

## Repo structure
```text
.
├─ fastq_filter.py             # Entry: imports + 2 facades
└─ modules/
   ├─ dna_rna_tools_module.py  # DNA/RNA helpers (pure functions)
   ├─ fastq_utils.py           # FASTQ filtering logic (GC%, length, phred33)
   └─ example_data.py          # (optional) EXAMPLE_FASTQ for local testing
```
## Install / dev setup

```bash
git clone https://github.com/alinatgrv/fastq-filter.git
cd fastq-filter
# optional:
# python3 -m venv .venv && source .venv/bin/activate
```

## Quick start

### DNA/RNA helpers (facade)
```python
from fastq_filter import run_dna_rna_tools

run_dna_rna_tools("ATGC", "reverse_complement")   # "GCAT"
run_dna_rna_tools("AUGC", "transcribe")           # "ATGC" (RNA → DNA)
run_dna_rna_tools("ATGC", "is_nucleic_acid")      # True
run_dna_rna_tools("ATGC", "GGG", "reverse")       # ["CGTA", "GGG"]
```

Supported: `is_nucleic_acid` | `transcribe` | `reverse` | `complement` | `reverse_complement`.

### FASTQ filtering
```python
from fastq_filter import filter_fastq

reads = {
    "@r1": ("ATGCATGC", "IIIIIIII"),  # len=8, GC=50%, meanQ≈40 (Phred+33)
    "@r2": ("AAAAAA",   "!!!!!!"),    # meanQ≈0
}
res = filter_fastq(
    reads,
    gc_bounds=(40, 60),      # inclusive
    length_bounds=(6, 10),   # inclusive
    quality_threshold=20.0   # mean Phred+33 >= 20
)
print(res.keys())
```

### With sample data:
```python
from modules.example_data import EXAMPLE_FASTQ
res = filter_fastq(EXAMPLE_FASTQ, gc_bounds=(40, 90), length_bounds=(60, 120), quality_threshold=20)
print(len(res))
```

## API

- `run_dna_rna_tools(*args: str) -> str | bool | list`  
  Last argument is the **procedure name**.  
  Returns a single value for one sequence, or a list for multiple sequences.

- `filter_fastq(seqs, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0.0) -> dict`  
  Parameters:
  - `seqs`: dict `name -> (sequence, quality_string)` (sequence and quality **lengths must match**).
  - `gc_bounds`: **inclusive** GC% interval. A single number means **upper bound** → `(0, number)`.
  - `length_bounds`: **inclusive** length interval. A single number means **upper bound** → `(0, number)`.
  - `quality_threshold`: inclusive lower bound on **mean** Phred+33 quality.  
    Phred+33 reminder: per-char `Q = ord(char) - 33`.

## References

- Illumina Quality Score Encoding (Phred+33): https://support.illumina.com/help/BaseSpace_Sequence_Hub_OLH_009008_2/Source/Informatics/BS/QualityScoreEncoding_swBS.htm  
- Python `ord()`: https://docs-python.ru/tutorial/vstroennye-funktsii-interpretatora-python/funktsija-ord/  
- ASCII table: https://www.asciitable.com/

## Development notes

- English docstrings; descriptive names (no one-letter variables; no leading “_” in function names).
- Modules keep pure logic; the entry file contains only imports and the two facades (`run_dna_rna_tools`, `filter_fastq`).
- All intervals are **inclusive**. Type hints throughout.

