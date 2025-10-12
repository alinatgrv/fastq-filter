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

def _has_mixed_tu(symbols: set) -> bool:
    return (("T" in symbols) or ("t" in symbols)) and (("U" in symbols) or ("u" in symbols))

def is_dna(seq: str) -> bool:
    """True if the sequence is valid DNA (ATCG only; no U; non-empty)."""
    if not seq:
        return False
    symbols = set(seq)
    if not symbols.issubset(VALID_DNA | VALID_RNA):
        return False
    if _has_mixed_tu(symbols):
        return False
    return symbols.issubset(VALID_DNA)

def is_rna(seq: str) -> bool:
    """True if the sequence is valid RNA (AUCG only; no T; non-empty)."""
    if not seq:
        return False
    symbols = set(seq)
    if not symbols.issubset(VALID_DNA | VALID_RNA):
        return False
    if _has_mixed_tu(symbols):
        return False
    return symbols.issubset(VALID_RNA)

def is_nucleic_acid(seq: str) -> bool:
    """True if the sequence is valid DNA or RNA (no T/U mixing)."""
    return is_dna(seq) or is_rna(seq)

def transcribe(seq: str) -> str:
    """
    Transcription preserving case:
    - DNA → RNA (T→U)
    - RNA → DNA (U→T)
    Raises ValueError for invalid sequence.
    """
    if is_dna(seq):
        return seq.replace("T", "U").replace("t", "u")
    if is_rna(seq):
        return seq.replace("U", "T").replace("u", "t")
    raise ValueError("Invalid sequence for transcription.")

def reverse(seq: str) -> str:
    """Return the reversed sequence. Validates input first."""
    if not is_nucleic_acid(seq):
        raise ValueError("Invalid nucleic acid sequence.")
    return seq[::-1]

def complement(seq: str) -> str:
    """Return the complementary sequence (DNA or RNA based on input)."""
    if is_dna(seq):
        table = COMPLEMENT_DNA
    elif is_rna(seq):
        table = COMPLEMENT_RNA
    else:
        raise ValueError("Invalid nucleic acid sequence.")
    return "".join(table.get(ch, ch) for ch in seq)

def reverse_complement(seq: str) -> str:
    """Return the reverse-complement sequence (DNA/RNA)."""
    return complement(reverse(seq))
