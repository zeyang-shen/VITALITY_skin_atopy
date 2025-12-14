import os
import shutil
import tempfile
import pytest
from unittest.mock import patch, call
from mgwas import core

@pytest.fixture
def temp_dir():
    path = tempfile.mkdtemp()
    yield path
    shutil.rmtree(path)

@pytest.fixture
def toy_data_dir():
    return os.path.join(os.path.dirname(__file__), "toy_data")

def test_parse_input_single(temp_dir):
    test_file = os.path.join(temp_dir, 'input.txt')
    with open(test_file, 'w') as f:
        f.write("sample\tphenotype\n")
        f.write("/path/to/sample1.fastq\t1\n")
        f.write("/path/to/sample2.fastq\t0\n")

    R1_list, R2_list, phenotype_dct = core.parse_input(test_file, paired=False)

    assert R1_list == ["/path/to/sample1.fastq", "/path/to/sample2.fastq"]
    assert R2_list == []
    assert "phenotype" in phenotype_dct
    assert phenotype_dct["phenotype"]["sample1.fastq"] == "1"
    assert phenotype_dct["phenotype"]["sample2.fastq"] == "0"


def test_parse_input_paired(temp_dir):
    test_file = os.path.join(temp_dir, 'input.txt')
    with open(test_file, 'w') as f:
        f.write("R1\tR2\tphenotype\n")
        f.write("/path/to/sample1_R1.fastq\t/path/to/sample1_R2.fastq\t1\n")
        f.write("/path/to/sample2_R1.fastq\t/path/to/sample2_R2.fastq\t0\n")

    R1_list, R2_list, phenotype_dct = core.parse_input(test_file, paired=True)

    assert R1_list == ["/path/to/sample1_R1.fastq", "/path/to/sample2_R1.fastq"]
    assert R2_list == ["/path/to/sample1_R2.fastq", "/path/to/sample2_R2.fastq"]
    assert "phenotype" in phenotype_dct
    assert phenotype_dct["phenotype"]["sample1_R1.fastq"] == "1"
    assert phenotype_dct["phenotype"]["sample2_R1.fastq"] == "0"


def test_prepare_reference_no_gbff(temp_dir, toy_data_dir):
    genome = os.path.join(toy_data_dir, 'reference.fna')

    ref_seq, gene_lst, bt2_index, snpeff_config = core.prepare_reference(
        genome=genome,
        genbank=None,
        temp=temp_dir,
        threads=1,
    )

    assert not snpeff_config
    assert os.path.exists(ref_seq)
    assert gene_lst == []
    assert os.path.exists(bt2_index+'.1.bt2')


def test_prepare_reference_with_snpeff(temp_dir, toy_data_dir):
    genome = os.path.join(toy_data_dir, 'reference.fna')
    gbff = os.path.join(toy_data_dir, 'reference.gbff')

    ref_seq, gene_lst, bt2_index, snpeff_config = core.prepare_reference(
        genome=genome,
        genbank=gbff,
        temp=temp_dir,
        threads=1,
    )
    
    assert os.path.exists(snpeff_config)
    assert os.path.exists(temp_dir+'/snpEff/data/reference/snpEffectPredictor.bin')
    assert os.path.exists(ref_seq)
    assert gene_lst[0][-1] == "EQW00_RS00005"
    assert os.path.exists(bt2_index+'.1.bt2')


def test_build_snpeff_and_annotate_variants(temp_dir, toy_data_dir):
    genome = os.path.join(toy_data_dir, 'reference.fna')
    gbff = os.path.join(toy_data_dir, 'reference.gbff')

    snpeff_config = core._build_snpeff_db(
        temp=temp_dir,
        seq_file=genome,
        gbk_file=gbff
    )
    
    assert os.path.exists(snpeff_config)
    assert os.path.exists(temp_dir+'/snpEff/data/reference/snpEffectPredictor.bin')
    
    vcf = os.path.join(toy_data_dir, 'output/all/combined.g.vcf')
    output_dir = os.path.join(temp_dir, "output")
    os.makedirs(output_dir, exist_ok=True)

    new_vcf, annotation_flag = core.annotate_variants(
        vcf=vcf,
        output_dir=output_dir,
        snpeff_config=snpeff_config
    )

    assert new_vcf == os.path.join(output_dir, "combined.g.annotated.vcf")
    assert os.path.exists(new_vcf)
    assert annotation_flag == True


# def test_prepare_reference_ncbi(temp_dir):
#     ref_seq, gene_lst, bt2_index = core.prepare_reference(
#         genome="Staphylococcus epidermidis",
#         genbank=None,
#         temp=temp_dir,
#         threads=1,
#     )

#     assert os.path.exists(ref_seq)
#     assert isinstance(gene_lst, list)
#     assert gene_lst[0][-1] == "EQW00_RS00005"
#     assert os.path.exists(bt2_index+'.1.bt2')


def test_map_reads_paired(temp_dir):
    R1_list = ["/path/to/sample1_R1.fastq.gz", "/path/to/sample2_R1.fastq.gz"]
    R2_list = ["/path/to/sample1_R2.fastq.gz", "/path/to/sample2_R2.fastq.gz"]
    paired = True
    bt2_index = "dummy_index"
    threads = 4

    with patch("mgwas.utils.file_exists", return_value=False) as mock_exists, \
         patch("mgwas.utils.run_command") as mock_run:

        core.map_reads(R1_list, R2_list, paired, bt2_index, temp_dir, threads)

        assert mock_exists.call_count == 2
        assert mock_run.call_count == 2

        # Verify the exact commands
        expected_calls = [
            call(f"bowtie2 --sensitive-local -x {bt2_index} -1 {R1_list[0]} -2 {R2_list[0]} -S {temp_dir}/sample1_R1.fastq.gz.sam --no-unal --threads 4 2>{temp_dir}/sample1_R1.fastq.gz.log"),
            call(f"bowtie2 --sensitive-local -x {bt2_index} -1 {R1_list[1]} -2 {R2_list[1]} -S {temp_dir}/sample2_R1.fastq.gz.sam --no-unal --threads 4 2>{temp_dir}/sample2_R1.fastq.gz.log")
        ]
        mock_run.assert_has_calls(expected_calls, any_order=False)


def test_map_reads_single(temp_dir):
    R1_list = ["/path/to/sample1_R1.fastq.gz", "/path/to/sample2_R1.fastq.gz"]
    R2_list = []  # not used
    paired = False
    bt2_index = "dummy_index"
    threads = 4

    with patch("mgwas.utils.file_exists", return_value=False) as mock_exists, \
         patch("mgwas.utils.run_command") as mock_run:

        core.map_reads(R1_list, R2_list, paired, bt2_index, temp_dir, threads)

        assert mock_exists.call_count == 2
        assert mock_run.call_count == 2

        # Verify the exact commands
        expected_calls = [
            call(f"bowtie2 --sensitive-local -x {bt2_index} -U {R1_list[0]} -S {temp_dir}/sample1_R1.fastq.gz.sam --no-unal --threads 4 2>{temp_dir}/sample1_R1.fastq.gz.log"),
            call(f"bowtie2 --sensitive-local -x {bt2_index} -U {R1_list[1]} -S {temp_dir}/sample2_R1.fastq.gz.sam --no-unal --threads 4 2>{temp_dir}/sample2_R1.fastq.gz.log")
        ]
        mock_run.assert_has_calls(expected_calls, any_order=False)


def test_process_coverage(temp_dir):
    threads = 4
    fake_samfiles = ["sample1.sam"]

    with patch("mgwas.utils.list_files_with_extension", return_value=fake_samfiles), \
         patch("mgwas.utils.file_exists", return_value=False), \
         patch("mgwas.utils.run_command") as mock_run:

        core.process_coverage(temp_dir, threads)

        prefix = os.path.join(temp_dir, "sample1")
        sam_path = f"{prefix}.sam"
        sorted_bam = f"{prefix}.sorted.bam"
        group_bam = f"{prefix}.sorted.groupAdded.bam"
        coverage = f"{prefix}.coverage"

        expected_calls = [
            call(f"samtools view -bS -@ {threads} {sam_path} | samtools sort -T {temp_dir} -@ {threads} > {sorted_bam}"),
            call(f"gatk AddOrReplaceReadGroups -I {sorted_bam} -O {group_bam} --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM sample1"),
            call(f"samtools index -@ {threads} {group_bam}"),
            call(f"samtools coverage {group_bam} > {coverage}"),
        ]

        mock_run.assert_has_calls(expected_calls, any_order=False)
        assert mock_run.call_count == 4


def test_analyze_sequencing_depths_binary(tmp_path):
    # Simulate a binary phenotype with 4 samples per group
    samples_high = ["sample1", "sample2", "sample3", "sample4"]
    samples_low  = ["sample5", "sample6", "sample7", "sample8"]

    phenotype_dct = {
        "pheno1": {
            **{s: "1" for s in samples_high},
            **{s: "0" for s in samples_low},
        }
    }

    # Coverage values for each sample
    coverage_values = {
        "sample1": 1000,
        "sample2": 1100,
        "sample3": 1100,
        "sample4": 800,
        "sample5": 100,
        "sample6": 90,
        "sample7": 80,
        "sample8": 100,
    }

    # Write dummy coverage files
    for sample, depth in coverage_values.items():
        file = tmp_path / f"{sample}.coverage"
        file.write_text(
            "# Example header\n"
            f"contig1\t100\t200\t{depth}\n"
        )

    # Run the function
    output_stats_file = tmp_path / "subsampling.stats"
    subsample_pheno_lst = core.analyze_sequencing_depths(
        temp_dir=tmp_path,
        phenotype_dct=phenotype_dct,
        subsample_stats_file=output_stats_file
    )

    # Assertions
    assert "pheno1" in subsample_pheno_lst
    assert output_stats_file.exists()

    lines = output_stats_file.read_text().strip().split("\n")
    samples_in_file = [line.split("\t")[1] for line in lines]

    # All high-depth samples should be listed for subsampling
    for sample in samples_high:
        assert sample in samples_in_file

    # Low-depth samples should not appear in the subsampling file
    for sample in samples_low:
        assert sample not in samples_in_file


def test_subsample_bam_files(tmp_path):
    output_dir = tmp_path / "output"
    phenotype_dct = {
        "phenotype1": {
            "sampleA": "1",
            "sampleB": "0",
            "sampleC": "1",
            "sampleD": "0"
        }
    }

    # Create a dummy `subsampling.stats` file
    subsampling_stats = tmp_path / "subsampling.stats"
    subsampling_stats.write_text(
        "phenotype1\tsampleA\t100\t1\t50\t0.5\n"
        "phenotype1\tsampleB\t80\t0\t40\t0.0\n"  # This should be skipped
        "phenotype1\tsampleC\t90\t1\t45\t0.01\n"
        "phenotype1\tsampleD\t90\t1\t45\t1.1\n"  # This should use original BAM
    )

    subsample_pheno_lst = ['phenotype1']

    # Create dummy BAM files
    for sample in ["sampleA", "sampleB", "sampleC", "sampleD"]:
        bam_file = tmp_path / f"{sample}.sorted.groupAdded.bam"
        bam_file.touch()

    # Mock run_command to prevent actual subprocess call
    with patch("mgwas.utils.run_command") as mock_run_cmd:
        bam_list_files = core.subsample_bam_files(
            temp_dir=tmp_path,
            output_dir=output_dir,
            phenotype_dct=phenotype_dct,
            subsample_pheno_lst=subsample_pheno_lst,
            subsample_stats_file=subsampling_stats,
            threads=4,
            subsample_seed=12
        )

        expected_calls = [
            call(
                f"samtools view -bs 12.5 {tmp_path / 'sampleA.sorted.groupAdded.bam'} > {tmp_path / 'sampleA.sorted.groupAdded.subsampled0.bam'}"
            ),
            call(
                f"samtools index -@ 4 {tmp_path / 'sampleA.sorted.groupAdded.subsampled0.bam'}"
            ),
            call(
                f"samtools view -bs 12.01 {tmp_path / 'sampleC.sorted.groupAdded.bam'} > {tmp_path / 'sampleC.sorted.groupAdded.subsampled0.bam'}"
            ),
            call(
                f"samtools index -@ 4 {tmp_path / 'sampleC.sorted.groupAdded.subsampled0.bam'}"
            ),
        ]
        mock_run_cmd.assert_has_calls(expected_calls, any_order=False)

    # Check output bam list file exists
    assert bam_list_files == [str(tmp_path / "bam.0.list")]
    
    bam_list_file = tmp_path / "bam.0.list"
    assert bam_list_file.exists()

    lines = bam_list_file.read_text().splitlines()
    assert str(tmp_path / "sampleA.sorted.groupAdded.subsampled0.bam") in lines
    assert str(tmp_path / "sampleC.sorted.groupAdded.subsampled0.bam") in lines
    assert str(tmp_path / "sampleD.sorted.groupAdded.bam") in lines  # should be added directly

