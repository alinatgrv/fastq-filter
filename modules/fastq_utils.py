from typing import Dict, Tuple, Union

FastqDict = Dict[str, Tuple[str, str]]
NumericOrBounds = Union[float, int, Tuple[Union[float, int], Union[float, int]]]

PHRED33_OFFSET = 33


def ensure_bounds_tuple(bounds: NumericOrBounds) -> Tuple[float, float]:
    """
    Normalize a bounds value to an inclusive (lower, upper) tuple.

    Rules
    -----
    - If `bounds` is a single number -> interpret as (0, bounds).
    - If `bounds` is a 2-tuple -> use as-is.
    - Order is normalized so that lower <= upper.

    Raises
    ------
    TypeError
        If bounds is neither a number nor a length-2 tuple of numbers.
    """
    if isinstance(bounds, (int, float)):
        lower, upper = 0.0, float(bounds)
    elif (
        isinstance(bounds, tuple)
        and len(bounds) == 2
        and all(isinstance(x, (int, float)) for x in bounds)
    ):
        lower, upper = float(bounds[0]), float(bounds[1])
    else:
        raise TypeError("Bounds must be a number or a (lower, upper) tuple of numbers.")
    if lower > upper:
        lower, upper = upper, lower
    return lower, upper


def compute_gc_percent(sequence: str) -> float:
    """Compute GC content (%) for a given sequence."""
    if not sequence:
        return 0.0
    seq = sequence.upper()
    gc_count = seq.count("G") + seq.count("C")
    return gc_count / len(seq) * 100


def compute_mean_phred33(quality_string: str) -> float:
    """
    Compute mean Phred+33 quality from a FASTQ quality string.

    Each character encodes a quality score Q via: Q = ord(char) - 33
    """
    if not quality_string:
        return 0.0
    total_quality = sum(ord(ch) - PHRED33_OFFSET for ch in quality_string)
    return total_quality / len(quality_string)


def value_in_inclusive_bounds(value: float, bounds: Tuple[float, float]) -> bool:
    """Return True if value is within inclusive bounds [lower, upper]."""
    lower, upper = bounds
    return lower <= value <= upper


def filter_reads(
    seqs: FastqDict,
    gc_bounds: NumericOrBounds = (0, 100),
    length_bounds: NumericOrBounds = (0, 2 ** 32),
    quality_threshold: Union[int, float] = 0
) -> FastqDict:
    """
    Filter FASTQ-like records by GC%, read length, and mean quality (Phred+33).
    """
    # Normalize bounds once
    gc_interval = ensure_bounds_tuple(gc_bounds)
    length_interval = ensure_bounds_tuple(length_bounds)
    q_threshold = float(quality_threshold)

    # Local predicates (closed over normalized bounds / threshold)
    def is_length_ok(sequence: str) -> bool:
        return value_in_inclusive_bounds(len(sequence), length_interval)

    def is_gc_ok(sequence: str) -> bool:
        return value_in_inclusive_bounds(compute_gc_percent(sequence), gc_interval)

    def is_quality_ok(quality: str) -> bool:
        return compute_mean_phred33(quality) >= q_threshold

    filtered_seqs: FastqDict = {}

    for name, pair in seqs.items():
        if not isinstance(pair, tuple) or len(pair) != 2:
            raise ValueError(f"Invalid record for '{name}': expected (sequence, quality_string)")

        sequence, quality = pair

        # Strict length consistency (sequence and quality must align)
        if len(sequence) != len(quality):
            raise ValueError(
                f"Length mismatch for '{name}': seq={len(sequence)}, qual={len(quality)}"
            )

        if is_length_ok(sequence) and is_gc_ok(sequence) and is_quality_ok(quality):
            filtered_seqs[name] = (sequence, quality)

    return filtered_seqs
