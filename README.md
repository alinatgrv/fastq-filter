# Bio Utilities: Sequence Classes & FASTQ Filter

A small collection of bioinformatics utilities implemented in Python.

The repository provides:

- Object-oriented biological sequence classes (DNA, RNA, proteins)
- FASTQ filtering using Biopython

---

## Installation

```bash
pip install -r requirements.txt
```

## Dependencies

- `biopython`

---

## Sequence Classes

Implemented in `fastq_filter.py`.

### Supported Types

- `DNASequence`
- `RNASequence`
- `AminoAcidSequence`

### Features

- Length support: `len(seq)`
- Indexing and slicing
- Alphabet validation
- Reverse / complement / reverse-complement (DNA/RNA)
- DNA → RNA transcription
- Amino acid composition calculation

### Example

```python
from fastq_filter import DNASequence

dna = DNASequence("ATGC")

print(dna.reverse_complement())
print(dna.transcribe())
```

---

## Sequence Classes

Implemented in `fastq_filter.py`.

### Supported Types

- `DNASequence`
- `RNASequence`
- `AminoAcidSequence`

### Features

- Length support: `len(seq)`
- Indexing and slicing
- Alphabet validation
- Reverse / complement / reverse-complement (DNA/RNA)
- DNA → RNA transcription
- Amino acid composition calculation

### Example

```python
from fastq_filter import DNASequence

dna = DNASequence("ATGC")

print(dna.reverse_complement())
print(dna.transcribe())
```

---

## FASTQ Filtering

### FASTQ files are filtered using Biopython:

- `SeqIO`

- `SeqRecord`

- `SeqUtils.gc_fraction`

### Filtering criteria

- GC percentage

- Read length

- Mean Phred+33 quality

### Usage

```python
python fastq_filter.py \
  --input_fastq example.fastq \
  --output_fastq filtered.fastq \
  --gc_bounds 30 70 \
  --length_bounds 50 300 \
  --quality_threshold 20 \
  --force
```

Output files are written to the `filtered/` directory.


---

### Notes

- Supports `.fastq` and `.fastq.gz`

- Output overwrite protection unless `--force` is specified
