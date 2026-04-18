import os
import logging
import pytest

from fastq_filter import (
    BiologicalSequence,
    DNASequence,
    RNASequence,
    AminoAcidSequence,
    open_maybe_gzip,
    filter_fastq,
    setup_logger,
)


@pytest.fixture
def sample_fastq(tmp_path):
    fastq_content = (
        "@read1\n"
        "ATGCATGC\n"
        "+\n"
        "IIIIIIII\n"
        "@read2\n"
        "AAAAAA\n"
        "+\n"
        "!!!!!!\n"
    )
    input_path = tmp_path / "sample.fastq"
    input_path.write_text(fastq_content)
    return input_path


class TestSequences:
    def test_biological_sequence_type_error(self):
        with pytest.raises(TypeError):
            BiologicalSequence(123)

    def test_dna_transcribe(self):
        dna = DNASequence("ATGC")
        assert dna.transcribe().sequence == "AUGC"

    def test_rna_invalid_symbol(self):
        with pytest.raises(ValueError):
            RNASequence("AUGTX")

    def test_amino_acid_composition(self):
        aa = AminoAcidSequence("AACC")
        composition = aa.aa_composition()
        assert composition["A"] == 0.5
        assert composition["C"] == 0.5


class TestFileIO:
    def test_open_maybe_gzip_write_and_read_plain_file(self, tmp_path):
        file_path = tmp_path / "test.txt"

        with open_maybe_gzip(file_path, "wt") as handle:
            handle.write("hello")

        with open_maybe_gzip(file_path, "rt") as handle:
            content = handle.read()

        assert content == "hello"


class TestFastqFilter:
    def test_filter_fastq_creates_output_file(self, sample_fastq, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        kept, total, output_path = filter_fastq(
            input_fastq_path=str(sample_fastq),
            output_fastq_name="filtered.fastq",
            quality_threshold=20,
        )

        assert total == 2
        assert kept == 1
        assert os.path.exists(output_path)

    def test_filter_fastq_filters_low_quality_reads(self, sample_fastq, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)

        kept, total, output_path = filter_fastq(
            input_fastq_path=str(sample_fastq),
            output_fastq_name="filtered.fastq",
            quality_threshold=20,
        )

        with open(output_path, "rt") as handle:
            content = handle.read()

        assert total == 2
        assert kept == 1
        assert "@read1" in content
        assert "@read2" not in content

    def test_filter_fastq_raises_if_output_exists(self, sample_fastq, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("filtered", exist_ok=True)

        existing_file = os.path.join("filtered", "filtered.fastq")
        with open(existing_file, "w") as handle:
            handle.write("existing")

        with pytest.raises(FileExistsError):
            filter_fastq(
                input_fastq_path=str(sample_fastq),
                output_fastq_name="filtered.fastq",
                quality_threshold=20,
            )

    def test_logging_writes_error_message(self, sample_fastq, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        os.makedirs("filtered", exist_ok=True)

        existing_file = os.path.join("filtered", "filtered.fastq")
        with open(existing_file, "w") as handle:
            handle.write("existing")

        logger = setup_logger("test.log")

        with pytest.raises(FileExistsError):
            filter_fastq(
                input_fastq_path=str(sample_fastq),
                output_fastq_name="filtered.fastq",
                quality_threshold=20,
                logger=logger,
            )

        with open("test.log", "r") as log_handle:
            log_content = log_handle.read()

        assert "ERROR" in log_content
        assert "already exists" in log_content