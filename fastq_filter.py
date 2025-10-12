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

# FASTQ helper utilities (no high-level filter function import)
from modules.fastq_utils import (
    ensure_bounds_tuple,
    compute_gc_percent,
    compute_mean_phred33,
    value_in_inclusive_bounds,
)

# ===== DNA/RNA processing =====

def run_dna_rna_tools(*args: str) -> Union[str, bool, List[Union[str, bool]]]:
    """
    Accepts one or more sequences (str) followed by a procedure name (str).

    Supported procedures:
      - 'is_nucleic_acid'     -> bool
      - 'transcribe'          -> str (DNA↔RNA)
      - 'reverse'             -> str
      - 'complement'          -> str
      - 'reverse_complement'  -> str

    Returns
    -------
    - bool or str for a single sequence
    - List[bool|str] for multiple sequences
    """
    if not args or len(args) < 2:
        raise ValueError("Provide at least one sequence and a procedure name.")

    *sequences, proc_name = args
    procedure = proc_name.strip().lower()

    # Function dispatch table instead of if-elif chain
    dispatch = {
        "is_nucleic_acid": is_nucleic_acid,
        "transcribe": transcribe,
        "reverse": reverse,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    func = dispatch.get(procedure)
    if func is None:
        supported = ", ".join(sorted(dispatch.keys()))
        raise ValueError(f"Unknown procedure: {proc_name!r}. Supported: {supported}")

    if len(sequences) == 1:
        return func(sequences[0])
    return [func(s) for s in sequences]


# ===== FASTQ processing =====

def filter_fastq(
    seqs: dict[str, tuple[str, str]],
    gc_bounds: Union[float, tuple[float, float]] = (0.0, 100.0),
    length_bounds: Union[int, tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0.0,
) -> dict[str, tuple[str, str]]:
    """
    Filter FASTQ reads by GC%, length, and mean quality (Phred+33), using inclusive bounds.

    This is the main high-level function that:
    - normalizes bounds,
    - validates read tuples (len(seq) == len(qual)),
    - computes metrics using helper utilities,
    - applies selection criteria,
    - collects and returns passing reads.
    """
    # 1) Normalize intervals (inclusive bounds)
    gc_interval = ensure_bounds_tuple(gc_bounds)          # -> (lower_gc, upper_gc)
    length_interval = ensure_bounds_tuple(length_bounds)  # -> (lower_len, upper_len)

    filtered: dict[str, tuple[str, str]] = {}

    # 2) Main loop over reads
    for read_name, pair in seqs.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError(f"Invalid record for '{read_name}': expected (sequence, quality_string)")

        sequence, quality_string = pair

        # 3) Validate length consistency
        if len(sequence) != len(quality_string):
            raise ValueError(
                f"Length mismatch for read '{read_name}': "
                f"sequence={len(sequence)}, quality={len(quality_string)}"
            )

        # 4) Compute metrics
        read_length = len(sequence)
        gc_percent = compute_gc_percent(sequence)                 # 0..100
        mean_quality = compute_mean_phred33(quality_string)       # average Q (Phred+33)

        # 5) Check criteria
        passes_gc = value_in_inclusive_bounds(gc_percent, gc_interval)
        passes_length = value_in_inclusive_bounds(read_length, length_interval)
        passes_quality = mean_quality >= float(quality_threshold)

        # 6) Keep read if it passes all filters
        if passes_gc and passes_length and passes_quality:
            filtered[read_name] = (sequence, quality_string)

    return filtered
