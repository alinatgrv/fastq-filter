"""
Entry point.

Contains ONLY:
- imports
- run_dna_rna_tools(...)
- filter_fastq(...)
"""

from typing import List, Union

# DNA/RNA helpers from the module
from modules.dna_rna_tools_module import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
)

# FASTQ filtering implementation
from modules.fastq_utils import filter_reads


# ===== DNA/RNA processing =====

ProcedureResult = Union[str, bool]
ManyResults = List[ProcedureResult]


def apply_procedure(proc: str, seq: str) -> ProcedureResult:
    """
    Apply a chosen DNA/RNA procedure to one sequence.

    Supported procedures:
      - 'is_nucleic_acid'     -> bool
      - 'transcribe'          -> str (DNA↔RNA)
      - 'reverse'             -> str
      - 'complement'          -> str
      - 'reverse_complement'  -> str
    """
    procedure = proc.strip().lower()
    if procedure == "is_nucleic_acid":
        return is_nucleic_acid(seq)
    if procedure == "transcribe":
        return transcribe(seq)
    if procedure == "reverse":
        return reverse(seq)
    if procedure == "complement":
        return complement(seq)
    if procedure == "reverse_complement":
        return reverse_complement(seq)
    raise ValueError(f"Unknown procedure: {proc!r}")


def run_dna_rna_tools(*args: str) -> Union[ProcedureResult, ManyResults]:
    """
    Accepts one or more sequences (str) followed by a procedure name (str).

    Examples
    --------
    run_dna_rna_tools("ATGC", "reverse_complement") -> "GCAT"
    run_dna_rna_tools("ATGC", "GGG", "reverse") -> ["CGTA", "GGG"]
    """
    if not args or len(args) < 2:
        raise ValueError("Provide at least one sequence and a procedure name.")
    *sequences, procedure = args
    if len(sequences) == 1:
        return apply_procedure(procedure, sequences[0])
    return [apply_procedure(procedure, s) for s in sequences]


# ===== FASTQ processing =====

def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: Union[float, tuple[float, float]] = (0.0, 100.0),
    length_bounds: Union[int, tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0.0,
) -> dict[str, tuple[str, str]]:
    """
    Filter FASTQ reads by GC%, length, and mean quality (phred33), inclusive.
    Delegates to modules.fastq_utils.filter_reads.
    """
    return filter_reads(
        seqs=seqs,
        gc_bounds=gc_bounds,
        length_bounds=length_bounds,
        quality_threshold=quality_threshold,
    )
