"""
Entry point.

Contains ONLY:
- imports
- run_dna_rna_tools(...)
- filter_fastq(...)
"""

from modules.in_out_utils import (
    iterate_fastq_records,
    build_safe_output_path,
    open_output_fastq,
    write_fastq_record,
)

from typing import List, Union

# DNA/RNA helpers from the module
from modules.dna_rna_tools_module import (
    is_nucleic_acid,
    transcribe,
    reverse,
    complement,
    reverse_complement,
)

# FASTQ helper utilities
from modules.fastq_utils import (
    compute_gc_percent,
    compute_mean_phred33,
    value_in_inclusive_bounds,
)

# ===== DNA/RNA processing =====

from typing import Iterable, Tuple

def normalize_bounds(bounds_value, default_lower, default_upper) -> Tuple[float, float]:
    """
    Normalize a scalar or a (lower, upper) sequence to a 2-tuple of floats.

    - If bounds_value is a single number -> treated as upper bound with lower=default_lower.
    - If bounds_value is a sequence of length >= 2 -> take the first two values.
    - Otherwise -> fall back to (default_lower, default_upper).
    """
    # sequence-like (list/tuple) with 1 or 2+ values
    if isinstance(bounds_value, (list, tuple)):
        if len(bounds_value) == 1:
            return float(default_lower), float(bounds_value[0])
        if len(bounds_value) >= 2:
            return float(bounds_value[0]), float(bounds_value[1])

    # scalar case
    try:
        upper_value = float(bounds_value)
        return float(default_lower), upper_value
    except (TypeError, ValueError):
        return float(default_lower), float(default_upper)


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
    input_fastq_path: str,
    output_fastq_name: str,
    gc_bounds: Union[float, tuple[float, float]] = (0.0, 100.0),
    length_bounds: Union[int, tuple[int, int]] = (0, 2**32),
    quality_threshold: float = 0.0,
    force_overwrite: bool = False,
) -> tuple[int, int, str]:
    """
    Stream-filter an input FASTQ and write passing records to an output FASTQ.

    The function:
    - normalizes bounds (GC% and length) using helper utilities,
    - computes metrics per record (GC%, length, mean Phred+33),
    - writes passing records immediately to avoid storing the whole file in memory.

    Parameters
    ----------
    input_fastq_path : str
        Path to the input FASTQ (or FASTQ.GZ) file.
    output_fastq_name : str
        File name for the output FASTQ written into the 'filtered' directory.
    gc_bounds : Union[float, tuple[float, float]], optional
        Inclusive GC% bounds. Either a single upper bound (lower=0) or a (lower, upper) tuple.
    length_bounds : Union[int, tuple[int, int]], optional
        Inclusive length bounds. Either a single upper bound (lower=0) or a (lower, upper) tuple.
    quality_threshold : float, optional
        Minimum mean Phred+33 quality to keep a read.
    force_overwrite : bool, optional
        If True, allow overwriting existing files. Otherwise a numeric suffix is added.

    Returns
    -------
    tuple[int, int, str]
        (kept_count, total_count, saved_output_path)
    """
    # 1) Normalize bounds using existing helpers from fastq_utils.py
    gc_interval = normalize_bounds(gc_bounds, 0.0, 100.0)
    length_interval = normalize_bounds(length_bounds, 0, 2**32)

    # 2) Prepare safe output path inside 'filtered/'
    output_path = build_safe_output_path(
        file_name=output_fastq_name,
        directory="filtered",
        force=force_overwrite,
    )

    kept_count = 0
    total_count = 0

    # 3) Stream records: read → check → write (if passed)
    with open_output_fastq(output_path) as output_handle:
        for header_line, sequence_line, plus_line, qualities_line in iterate_fastq_records(input_fastq_path):
            total_count += 1

            read_length = len(sequence_line)
            gc_percent = compute_gc_percent(sequence_line)
            mean_quality = compute_mean_phred33(qualities_line)

            passes_gc = value_in_inclusive_bounds(gc_percent, gc_interval)
            passes_length = value_in_inclusive_bounds(read_length, length_interval)
            passes_quality = mean_quality >= float(quality_threshold)

            if passes_gc and passes_length and passes_quality:
                write_fastq_record(
                    output_handle,
                    (header_line, sequence_line, plus_line, qualities_line),
                )
                kept_count += 1

    return kept_count, total_count, str(output_path)

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Filter FASTQ and write passing records to 'filtered/'."
    )
    parser.add_argument(
        "--input_fastq",
        required=True,
        help="Path to input FASTQ or FASTQ.GZ file.",
    )
    parser.add_argument(
        "--output_fastq",
        required=True,
        help="Output FASTQ file name 'filtered/'.",
    )
    parser.add_argument(
        "--gc_bounds",
        nargs="+",
        type=float,
        default=[0.0, 100.0],
        help="GC%% bounds: give one value (upper) or two values (lower upper).",
    )
    parser.add_argument(
        "--length_bounds",
        nargs="+",
        type=int,
        default=[0, 4294967296],
        help="Length bounds: one value (upper) or two values (lower upper).",
    )
    parser.add_argument(
        "--quality_threshold",
        type=float,
        default=0.0,
        help="Minimum mean Phred+33 quality to keep a read.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting existing output file name.",
    )

    args = parser.parse_args()

    if len(args.gc_bounds) == 1:
        gc_arg = args.gc_bounds[0]
    else:
        gc_arg = (args.gc_bounds[0], args.gc_bounds[1])

    if len(args.length_bounds) == 1:
        length_arg = args.length_bounds[0]
    else:
        length_arg = (args.length_bounds[0], args.length_bounds[1])

    kept_count, total_count, saved_path = filter_fastq(
        input_fastq_path=args.input_fastq,
        output_fastq_name=args.output_fastq,
        gc_bounds=gc_arg,
        length_bounds=length_arg,
        quality_threshold=args.quality_threshold,
        force_overwrite=args.force,
    )
    print(f"Kept {kept_count} / {total_count} reads. Saved to: {saved_path}")
