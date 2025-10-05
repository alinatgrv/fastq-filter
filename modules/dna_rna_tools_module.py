"""
dna_rna_tools_module.py
Utilities for working with DNA/RNA sequences.
"""

from typing import Union

# Allowed alphabets
VALID_DNA = set("ATCGatcg")
VALID_RNA = set("AUCGaucg")

# Complementarity tables
COMPLEMENT_DNA = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "a": "t", "t": "a", "c": "g", "g": "c",
}
COMPLEMENT_RNA = {
    "A": "U", "U": "A", "C": "G", "G": "C",
    "a": "u", "u": "a", "c": "g", "g": "c",
}


def detect_seq_type(seq: str) -> Union[str, None]:
    """
    Returns:
      - 'DNA' if the sequence contains only ATCG,
      - 'RNA' if the sequence contains only AUCG,
      - None if it has invalid symbols, mixes T and U, or is empty.
    """
    if not seq:
        return None

    symbols = set(seq)

    # Only DNA/RNA alphabet?
    if not symbols.issubset(VALID_DNA | VALID_RNA):
        return None

    # Do not allow T and U simultaneously
    has_t = "T" in symbols or "t" in symbols
    has_u = "U" in symbols or "u" in symbols
    if has_t and has_u:
        return None

    if symbols.issubset(VALID_DNA):
        return "DNA"
    if symbols.issubset(VALID_RNA):
        return "RNA"
    return None


def is_nucleic_acid(seq: str) -> bool:
    """True if the sequence is a valid DNA or RNA (no T/U mixing)."""
    return detect_seq_type(seq) in ("DNA", "RNA")


def transcribe(seq: str) -> str:
    """
    Transcription preserving case:
    - DNA → RNA (T→U)
    - RNA → DNA (U→T)
    Raises ValueError for invalid sequence.
    """
    seq_type = detect_seq_type(seq)
    if seq_type is None:
        raise ValueError("Invalid sequence for transcription.")
    if seq_type == "DNA":
        return seq.replace("T", "U").replace("t", "u")
    return seq.replace("U", "T").replace("u", "t")


def reverse(seq: str) -> str:
    """Return the reversed sequence. Validates input first."""
    if not is_nucleic_acid(seq):
        raise ValueError("Invalid nucleic acid sequence.")
    return seq[::-1]


def complement(seq: str) -> str:
    """Return the complementary sequence (DNA or RNA based on input)."""
    seq_type = detect_seq_type(seq)
    if seq_type is None:
        raise ValueError("Invalid nucleic acid sequence.")
    table = COMPLEMENT_DNA if seq_type == "DNA" else COMPLEMENT_RNA
    return "".join(table.get(ch, ch) for ch in seq)


def reverse_complement(seq: str) -> str:
    """Return the reverse-complement sequence."""
    seq_type = detect_seq_type(seq)
    if seq_type is None:
        raise ValueError("Invalid nucleic acid sequence.")
    table = COMPLEMENT_DNA if seq_type == "DNA" else COMPLEMENT_RNA
    return "".join(table.get(ch, ch) for ch in reversed(seq))
