from typing import Optional

# Допустимые символы для последовательностей
VALID_DNA = set("ATCGatcg")
VALID_RNA = set("AUCGaucg")

# Таблицы комплементарности
COMPLEMENT_DNA = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C",
    "a": "t",
    "t": "a",
    "c": "g",
    "g": "c",
}
COMPLEMENT_RNA = {
    "A": "U",
    "U": "A",
    "C": "G",
    "G": "C",
    "a": "u",
    "u": "a",
    "c": "g",
    "g": "c",
}

def _seq_type(seq: str) -> Optional[str]:
    """
    'DNA' если только ATCG,
    'RNA' если только AUCG,
    None если есть посторонние символы или одновременно T и U.
    """
    if not seq:  # строка пустая
        return None

    symbols = set(seq)
    # посторонние символы
    if not symbols.issubset(VALID_DNA | VALID_RNA):
        return None
    # смешение T и U недопустимо
    if (("T" in symbols or "t" in symbols) and ("U" in symbols or "u" in symbols)):
        return None
    if symbols.issubset(VALID_DNA):
        return "DNA"
    if symbols.issubset(VALID_RNA):
        return "RNA"
    return None

def is_nucleic_acid(seq: str) -> bool:
    """Возвращает True, если последовательность валидна (DNA или RNA)."""
    return _seq_type(seq) in ("DNA", "RNA")

def transcribe(seq: str) -> str:
    """
    DNA → RNA (T→U), RNA → DNA (U→T). Регистр сохраняется.
    Бросает ValueError при невалидной последовательности.
    """
    st = _seq_type(seq)
    if st is None:
        raise ValueError("Invalid sequence for transcription.")
    if st == "DNA":
        return seq.replace("T", "U").replace("t", "u")
    # RNA
    return seq.replace("U", "T").replace("u", "t")

def reverse(seq: str) -> str:
    """Возвращает развёрнутую последовательность (с проверкой валидности)."""
    if not is_nucleic_acid(seq):
        raise ValueError("Invalid nucleic acid sequence.")
    return seq[::-1]

def complement(seq: str) -> str:
    """Возвращает комплементарную последовательность."""
    st = _seq_type(seq)
    if st is None:
        raise ValueError("Invalid nucleic acid sequence.")
    table = COMPLEMENT_DNA if st == "DNA" else COMPLEMENT_RNA
    return "".join(table.get(ch, ch) for ch in seq)

def reverse_complement(seq: str) -> str:
    """Возвращает обратную комплементарную последовательность."""
    st = _seq_type(seq)
    if st is None:
        raise ValueError("Invalid nucleic acid sequence.")
    table = COMPLEMENT_DNA if st == "DNA" else COMPLEMENT_RNA
    return "".join(table.get(ch, ch) for ch in reversed(seq))
