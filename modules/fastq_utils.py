from typing import Dict, Tuple, Union

FastqDict = Dict[str, Tuple[str, str]]
NumericOrBounds = Union[float, int, Tuple[Union[float, int], Union[float, int]]]


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
    """
    Compute GC content percentage for a nucleotide sequence.

    Returns
    -------
    float
        GC% in [0, 100]. Empty sequence -> 0.0.
    """
    if not sequence:
        return 0.0
    uppercase_seq = sequence.upper()
    gc_count = sum(1 for base in uppercase_seq if base in ("G", "C"))
    return 100.0 * gc_count / len(uppercase_seq)


def compute_mean_phred33(quality_string: str) -> float:
    """
    Compute mean Phred+33 quality from a FASTQ quality string.

    Each character encodes a quality score Q via:  Q = ord(char) - 33
    """
    if not quality_string:
        return 0.0
    total_quality = sum(ord(ch) - 33 for ch in quality_string)
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

    Input format (FASTQ-like dict)
    ------------------------------
    {
      "read_name": ("SEQUENCE", "QUALITY_STRING"),
      ...
    }

    Filters (all conditions must be satisfied)
    -----------------------------------------
    1) GC% within `gc_bounds` (inclusive).
       - If a single number is given, it is treated as upper bound: (0, number).
    2) Read length within `length_bounds` (inclusive).
       - If a single number is given, it is treated as upper bound: (0, number).
    3) Mean quality (Phred+33) >= `quality_threshold`.

    Returns
    -------
    dict[name -> (sequence, quality)]

    Raises
    ------
    ValueError
        If a read has sequence and quality of different lengths.
    TypeError
        If bounds are provided in unsupported formats.
    """
    gc_interval = ensure_bounds_tuple(gc_bounds)
    length_interval = ensure_bounds_tuple(length_bounds)

    filtered: FastqDict = {}

    for read_name, (sequence, quality_string) in seqs.items():
        if len(sequence) != len(quality_string):
            raise ValueError(
                f"Length mismatch for read '{read_name}': "
                f"sequence={len(sequence)}, quality={len(quality_string)}"
            )

        read_length = len(sequence)
        gc_percent = compute_gc_percent(sequence)
        mean_quality = compute_mean_phred33(quality_string)

        passes_gc = value_in_inclusive_bounds(gc_percent, gc_interval)
        passes_length = value_in_inclusive_bounds(read_length, length_interval)
        passes_quality = mean_quality >= float(quality_threshold)

        if passes_gc and passes_length and passes_quality:
            filtered[read_name] = (sequence, quality_string)

    return filtered
