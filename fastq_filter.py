from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
import gzip
import argparse 

# 1: Sequences OOP

class BiologicalSequence:
    """
    Class for biological sequences.

    Provides:
    - len(sequence)
    - indexing and slicing
    - alphabet validation

    """

    def __init__(self, sequence):
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string.")

        self._seq = sequence
        self.check_alphabet()

    def __len__(self):
        """Return sequence length."""
        return len(self._seq)

    def __getitem__(self, item):
        """
        Support indexing and slicing.

        - If index -> return a single character (str)
        - If slice -> return a new object of the same class
        """
        if isinstance(item, slice):
            return self.__class__(self._seq[item])
        return self._seq[item]

    def __str__(self):
        """Readable string representation."""
        return f"{self.__class__.__name__}('{self._seq}')"

    def __repr__(self):
        """Official string representation."""
        return self.__str__()

    def check_alphabet(self):
        """
        Validate alphabet correctness.

        Must be implemented in subclasses.
        """
        raise NotImplementedError("Subclasses must implement check_alphabet().")

    @property
    def sequence(self):
        """Return raw sequence string."""
        return self._seq


class NucleicAcidSequence(BiologicalSequence):
    """
    Base class for nucleic acid sequences (DNA and RNA).

    Implements:
    - alphabet validation
    - complement()
    - reverse()
    - reverse_complement()

    Polymorphism realisation:
    Subclasses define their own alphabet and complement mapping.
    """

    _alphabet = None
    _complement_map = None

    def _ensure_configuration(self):
        """
        Ensure that subclass defined required attributes.

        Prevents direct usage of NucleicAcidSequence.
        """
        if self._alphabet is None or self._complement_map is None:
            raise NotImplementedError(
                "Subclasses must define _alphabet and _complement_map."
            )

    def check_alphabet(self):
        """Validate that sequence contains only allowed symbols."""
        self._ensure_configuration()

        for symbol in self._seq:
            if symbol not in self._alphabet:
                raise ValueError(
                    f"Invalid symbol '{symbol}' for {self.__class__.__name__}."
                )
        return True

    def complement(self):
        self._ensure_configuration()

        complemented = []
        for symbol in self._seq:
            complemented.append(self._complement_map[symbol])

        return self.__class__("".join(complemented))

    def reverse(self):
        return self.__class__(self._seq[::-1])

    def reverse_complement(self):
        return self.reverse().complement()


class DNASequence(NucleicAcidSequence):
    """
    DNA sequence class.
    """

    _alphabet = set("ATCGatcg")

    _complement_map = {
        "A": "T", "T": "A", "C": "G", "G": "C",
        "a": "t", "t": "a", "c": "g", "g": "c",
    }

    def transcribe(self):
        """
        Transcribe DNA to RNA (T -> U).
        """
        rna_sequence = self._seq.replace("T", "U").replace("t", "u")
        return RNASequence(rna_sequence)


class RNASequence(NucleicAcidSequence):
    """
    RNA sequence class.
    """

    _alphabet = set("AUCGaucg")

    _complement_map = {
        "A": "U", "U": "A", "C": "G", "G": "C",
        "a": "u", "u": "a", "c": "g", "g": "c",
    }


class AminoAcidSequence(BiologicalSequence):
    """
    Amino acid sequence class.

    Supports validation of standard 20 amino acids.
    """

    _alphabet = set("ACDEFGHIKLMNPQRSTVWY")

    def check_alphabet(self):
        """Validate amino acid alphabet."""
        for symbol in self._seq:
            if symbol not in self._alphabet:
                raise ValueError(
                    f"Invalid symbol '{symbol}' for AminoAcidSequence."
                )
        return True

    def aa_composition(self):
        """
        Calculate amino acid composition.

        Returns: 
        dict
            Dictionary of amino acid frequencies (fractions).
        """
        if len(self._seq) == 0:
            return {}

        counts = {}
        for symbol in self._seq:
            counts[symbol] = counts.get(symbol, 0) + 1

        total_length = len(self._seq)

        frequencies = {}
        for aa in counts:
            frequencies[aa] = counts[aa] / total_length

        return frequencies



# 2: FASTQ filter with Biopython


def open_maybe_gzip(path, mode):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8")
    return open(path, mode, encoding="utf-8")


def filter_fastq(
    input_fastq_path,
    output_fastq_name,
    gc_bounds=(0.0, 100.0),
    length_bounds=(0, 2**32),
    quality_threshold=0.0,
    force_overwrite=False,
):
    """
    Filter FASTQ using Biopython objects (SeqIO, SeqRecord, SeqUtils).
    """

    gc_lower, gc_upper = gc_bounds
    len_lower, len_upper = length_bounds

    os.makedirs("filtered", exist_ok=True)
    output_path = os.path.join("filtered", output_fastq_name)

    if os.path.exists(output_path) and not force_overwrite:
        raise FileExistsError(
            f"{output_path} already exists. Use --force to overwrite."
        )

    kept = 0
    total = 0

    with open_maybe_gzip(input_fastq_path, "rt") as in_handle, \
         open_maybe_gzip(output_path, "wt") as out_handle:

        for record in SeqIO.parse(in_handle, "fastq"):
            total += 1

            seq_length = len(record.seq)
            gc_percent = gc_fraction(record.seq) * 100
            qualities = record.letter_annotations["phred_quality"]
            mean_quality = sum(qualities) / len(qualities)

            if (
                gc_lower <= gc_percent <= gc_upper
                and len_lower <= seq_length <= len_upper
                and mean_quality >= quality_threshold
            ):
                SeqIO.write(record, out_handle, "fastq")
                kept += 1

    return kept, total, output_path

# CLI

def main():
    parser = argparse.ArgumentParser(
        description="Filter FASTQ/FASTQ.GZ and write passing reads to 'filtered/'."
    )
    parser.add_argument("--input_fastq", required=True, help="Input FASTQ or FASTQ.GZ path.")
    parser.add_argument("--output_fastq", required=True, help="Output file name (inside filtered/).")

    parser.add_argument(
        "--gc_bounds",
        nargs="+",
        type=float,
        default=[0.0, 100.0],
        help="GC%% bounds: provide 1 value (upper) or 2 values (lower upper).",
    )
    parser.add_argument(
        "--length_bounds",
        nargs="+",
        type=int,
        default=[0, 2**32],
        help="Length bounds: provide 1 value (upper) or 2 values (lower upper).",
    )
    parser.add_argument(
        "--quality_threshold",
        type=float,
        default=0.0,
        help="Minimum mean Phred+33 quality to keep a read.",
    )
    parser.add_argument("--force", action="store_true", help="Allow overwriting output file.")

    args = parser.parse_args()

    gc_arg = args.gc_bounds[0] if len(args.gc_bounds) == 1 else (args.gc_bounds[0], args.gc_bounds[1])
    length_arg = args.length_bounds[0] if len(args.length_bounds) == 1 else (args.length_bounds[0], args.length_bounds[1])

    kept, total, out_path = filter_fastq(
        input_fastq_path=args.input_fastq,
        output_fastq_name=args.output_fastq,
        gc_bounds=gc_arg,
        length_bounds=length_arg,
        quality_threshold=args.quality_threshold,
        force_overwrite=args.force,
    )

    print(f"Kept {kept} / {total} reads. Saved to: {out_path}")


if __name__ == "__main__":
    main()