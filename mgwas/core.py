import os
import logging
import shutil
from Bio import SeqIO
from collections import Counter, defaultdict
from scipy import stats
import numpy as np
import ncbi_genome_download as ngd
from mgwas import utils, gwas

def prepare_reference(genome, genbank, temp, threads):
    """Prepare reference genome, annotations, and indices."""
    ref_fna = os.path.join(temp, 'reference.fna')
    ref_gbff = os.path.join(temp, 'reference.gbff')
    gene_lst = []
    snpeff_config = None

    if not utils.file_exists(ref_fna):
        if utils.file_exists(genome):
            logging.info(f"Custom genome specified: {genome}")
            shutil.copy(genome, ref_fna)
            if genbank and utils.file_exists(genbank):
                logging.info(f"Gene annotation file specified: {genbank}")
                shutil.copy(genbank, ref_gbff)
            else:
                logging.info("No gene annotation file found")
        else:
            logging.info(f"Downloading reference genome '{genome}' from NCBI...")
            ngd.download(genera=genome, groups='bacteria', refseq_categories='reference', file_formats='fasta,genbank', output=temp)
            utils.move_to_root_folder(temp, temp)
            utils.gunzip_all(temp)
            utils.rename_file_given_extension(temp, 'fna')
            utils.rename_file_given_extension(temp, 'gbff')
    else:
        logging.info("Reference genome exists. Continue analysis...")

    # Load gene annotations only if file exists
    if utils.file_exists(ref_gbff):
        gene_lst = _load_gene_annotations(ref_gbff)
        snpeff_config = _build_snpeff_db(temp, ref_fna, ref_gbff)
    else:
        logging.info("Skipping annotation-related steps (No gene annotation file found).")

    # Prepare reference for GATK
    if not utils.file_exists(os.path.join(temp, 'reference.dict')):
        logging.info("Preparing reference for GATK...")
        utils.run_command(f"samtools faidx {ref_fna}")
        utils.run_command(f"gatk CreateSequenceDictionary -R {ref_fna}")

    # Create Bowtie2 index
    logging.info("Creating Bowtie2 index...")
    utils.run_command(f"bowtie2-build --threads {threads} {ref_fna} {temp}/reference")

    return ref_fna, gene_lst, temp+'/reference', snpeff_config


def _build_snpeff_db(temp, seq_file, gbk_file):
    logging.info("Creating snpEff database...")
    db_dir = os.path.join(temp, 'snpEff/data/reference/')
    os.makedirs(db_dir, exist_ok=True)

    # Copy files
    for fname, target in [(seq_file, 'sequences.fa'), (gbk_file, 'genes.gbk')]:
        target_path = os.path.join(db_dir, target)
        if utils.file_exists(target_path):
            os.remove(target_path)
        shutil.copy(fname, target_path)

    # Write configuration file
    config_path = os.path.join(temp, 'snpEff/snpEff.config')
    with open(config_path, 'w') as wf:
        wf.write(f'data.dir = ./data\n')
        wf.write(f'reference.genome: Reference Genome\n')

    # Create database
    utils.run_command(
        f"snpEff build -config {config_path} -genbank -v reference"
    )
    return config_path


def _load_gene_annotations(gbff_path):
    genes = []
    with open(gbff_path, "r") as handle:
        for record in SeqIO.parse(handle, "genbank"):
            seq_id = record.id
            for feature in record.features:
                if feature.type == "CDS":
                    try:
                        for part in feature.location.parts:
                            genes.append((seq_id, part.start, part.end, feature.qualifiers['locus_tag'][0]))
                    except:
                        continue
    return genes


def parse_input(input_file, paired):
    with open(input_file, 'r') as f:
        lines = f.readlines()

    if paired:
        phenotype_names = lines[0].strip().split('\t')[2:]
    else:
        phenotype_names = lines[0].strip().split('\t')[1:]

    phenotype_dct = {}
    R1_list, R2_list = [], []

    for line in lines[1:]:
        fields = line.strip().split('\t')
        sample = fields[0].split('/')[-1]
        if paired:
            R1_list.append(fields[0])
            R2_list.append(fields[1])
            phenos = fields[2:]
        else:
            R1_list.append(fields[0])
            phenos = fields[1:]

        for name, val in zip(phenotype_names, phenos):
            if val == '':
                continue
            phenotype_dct.setdefault(name, {})[sample] = val

    return R1_list, R2_list, phenotype_dct


def map_reads(R1_list, R2_list, paired, bt2_index, temp_dir, threads):
    """
    Map sequencing reads with bowtie2.
    """
    if paired:
        logging.info("Mapping paired-end data...")
        for R1, R2 in zip(R1_list, R2_list):
            prefix = os.path.basename(R1)
            sam_file = os.path.join(temp_dir, f"{prefix}.sam")
            log_file = os.path.join(temp_dir, f"{prefix}.log")
            if utils.file_exists(sam_file):
                logging.info(f"Skipping: {R1}. SAM file exists.")
                continue
            logging.info(f"Mapping: {R1}")
            try:
                utils.run_command(f"bowtie2 --sensitive-local -x {bt2_index} -1 {R1} -2 {R2} -S {sam_file} --no-unal --threads {threads} 2>{log_file}")
            except Exception:
                logging.info(f"{R1} failed mapping!")
    else:
        logging.info("Mapping single-read data...")
        for R1 in R1_list:
            prefix = os.path.basename(R1)
            sam_file = os.path.join(temp_dir, f"{prefix}.sam")
            log_file = os.path.join(temp_dir, f"{prefix}.log")
            if utils.file_exists(sam_file):
                logging.info(f"Skipping: {R1}. SAM file exists.")
                continue
            logging.info(f"Mapping: {R1}")
            try:
                utils.run_command(f"bowtie2 --sensitive-local -x {bt2_index} -U {R1} -S {sam_file} --no-unal --threads {threads} 2>{log_file}")
            except Exception:
                logging.info(f"{R1} failed mapping!")


def process_coverage(temp_dir, threads):
    """
    Convert SAM to BAM, add read groups, index, and compute coverage.
    """
    logging.info("Computing coverage...")
    samfiles = utils.list_files_with_extension(temp_dir, "sam")
    for samfile in samfiles:
        coverage_file = os.path.join(temp_dir, samfile.replace('.sam', '.coverage'))
        if utils.file_exists(coverage_file):
            logging.info(f"Skipping: {samfile}. Coverage file exists.")
            continue
        logging.info(f"Processing: {samfile}")
        
        # Convert SAM to BAM
        bamfile = samfile.replace('.sam', '.sorted.bam')
        bam_path = os.path.join(temp_dir, bamfile)
        sam_path = os.path.join(temp_dir, samfile)
        utils.run_command(f"samtools view -bS -@ {threads} {sam_path} | samtools sort -T {temp_dir} -@ {threads} > {bam_path}")
        
        # Convert to GATK format
        gatkfile = bamfile.replace('.bam', '.groupAdded.bam')
        gatk_path = os.path.join(temp_dir, gatkfile)
        utils.run_command(f"gatk AddOrReplaceReadGroups -I {bam_path} -O {gatk_path} --RGLB lib1 --RGPL ILLUMINA --RGPU unit1 --RGSM {samfile.replace('.sam', '')}")
        utils.run_command(f"samtools index -@ {threads} {gatk_path}")
        
        # Compute coverage
        utils.run_command(f"samtools coverage {gatk_path} > {coverage_file}")


def analyze_sequencing_depths(temp_dir, phenotype_dct, subsample_stats_file):
    """
    Analyze sequencing depths to determine if subsampling is needed.
    
    Parameters:
        temp_dir (str): Path to temporary directory with coverage files.
        phenotype_dct (dict): Mapping from phenotype names to sample:value dictionaries.
        subsample_stats_file (str): Output path to write subsampling stats.
    
    Returns:
        subsample_pheno_lst (list): List of phenotype names that need subsampling.
    """
    logging.info("Analyzing sequencing depths...")
    coverage_files = utils.list_files_with_extension(temp_dir, "coverage")
    reads_dct = {}
    for coverage_file in coverage_files:
        sample = coverage_file.replace('.coverage', '')
        total_reads = 0
        with open(os.path.join(temp_dir, coverage_file), 'r') as rf:
            for line in rf:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                total_reads += int(fields[3])
        reads_dct[sample] = total_reads

    subsample_pheno_lst = []
    with open(subsample_stats_file, 'w') as wf:
        for pheno, sample_val_dct in phenotype_dct.items():
            values = list(sample_val_dct.values())
            unique_vals = list(set(values))
            
            if len(unique_vals) > 2:
                # Quantitative phenotype
                samples = list(sample_val_dct.keys())
                reads = [reads_dct.get(s, 0) for s in samples]
                pheno_vals = [float(sample_val_dct[s]) for s in samples]
                corr, pval = stats.spearmanr(np.log2(np.array(reads)+1), pheno_vals)
                if pval < 0.05 and abs(corr) > 0.2:
                    logging.info(f"Phenotype {pheno} shows correlated read depth; subsampling recommended.")
                    subsample_pheno_lst.append(pheno)
                    min_read = min(reads)
                    for s, r, p in zip(samples, reads, pheno_vals):
                        sub_p = round(1 - ((r - min_read) / r), 4)
                        wf.write(f"{pheno}\t{s}\t{r}\t{p}\t{min_read}\t{sub_p}\n")

            elif len(unique_vals) == 2:
                # Binary phenotype
                group0, group1 = unique_vals
                group0_samples = [s for s, v in sample_val_dct.items() if v == group0]
                group1_samples = [s for s, v in sample_val_dct.items() if v == group1]
                r0 = np.array([reads_dct[s] for s in group0_samples])
                r1 = np.array([reads_dct[s] for s in group1_samples])
                try:
                    stat, pval = stats.mannwhitneyu(r0, r1)
                except ValueError:
                    continue
                if pval < 0.05:
                    if np.median(r0) > 2 * np.median(r1):
                        logging.info(f"Phenotype {pheno} shows imbalanced depths between groups: {group0} > {group1}")
                        subsample_pheno_lst.append(pheno)
                        for s, r in zip(group0_samples, r0):
                            p = np.sum(r0 < r) / len(r0) * 100
                            r_new = np.percentile(r1, p)
                            sub_p = round(1 - ((r - r_new) / r), 4)
                            wf.write(f"{pheno}\t{s}\t{r}\t{p}\t{r_new}\t{sub_p}\n")
                    elif np.median(r1) > 2 * np.median(r0):
                        logging.info(f"Phenotype {pheno} shows imbalanced depths between groups: {group1} > {group0}")
                        subsample_pheno_lst.append(pheno)
                        for s, r in zip(group1_samples, r1):
                            p = np.sum(r1 < r) / len(r1) * 100
                            r_new = np.percentile(r0, p)
                            sub_p = round(1 - ((r - r_new) / r), 4)
                            wf.write(f"{pheno}\t{s}\t{r}\t{p}\t{r_new}\t{sub_p}\n")

    return subsample_pheno_lst


def subsample_bam_files(
    temp_dir,
    output_dir,
    phenotype_dct,
    subsample_pheno_lst,
    subsample_stats_file,
    threads,
    subsample_seed=20
):
    """
    Subsample BAM files for phenotypes with imbalanced sequencing depth.
    
    Args:
        temp_dir (str): Temporary working directory.
        output_dir (str): Output directory to store phenotype maps.
        phenotype_dct (dict): Phenotype-to-sample mappings.
        subsample_pheno_lst (list): Phenotypes that require subsampling.
        subsample_stats_file (str): Path to 'subsampling.stats' file.
        threads (int): Number of threads for samtools.
        subsample_seed (int): Random seed prefix for samtools subsampling (default: 20).
    """
    logging.info("Performing BAM subsampling...")

    subsample_pheno_dct = {p: i for i, p in enumerate(subsample_pheno_lst)}

    # Write out phenotype index file
    os.makedirs(output_dir, exist_ok=True)
    with open(os.path.join(output_dir, "phenotype_subsampled.txt"), "w") as wf:
        for p, i in subsample_pheno_dct.items():
            wf.write(f"{p}\t{i}\n")

    # Track which samples need subsampling and their target files
    subsample_pheno_sample_dct = defaultdict(list)
    pheno_bam_file_dct = defaultdict(list)

    with open(subsample_stats_file, "r") as rf:
        for line in rf:
            pheno, sample, _, _, _, subsample_prop = line.strip().split('\t')
            subsample_prop = float(subsample_prop)

            # skip if subsample ratio is <= 0
            if subsample_prop <= 0:
                continue

            subsample_pheno_sample_dct[pheno].append(sample)
            
            if subsample_prop >= 1:
                pheno_bam_file_dct[pheno].append(f"{sample}.sorted.groupAdded.bam")
            else:
                subsample_str = str(subsample_prop)[1:]  # remove leading '0'
                new_bam = f"{sample}.sorted.groupAdded.subsampled{subsample_pheno_dct[pheno]}.bam"
                input_bam = os.path.join(temp_dir, f"{sample}.sorted.groupAdded.bam")
                output_bam = os.path.join(temp_dir, new_bam)
                logging.info(f"Subsampling {sample} by {subsample_str} for phenotype {pheno}...")
                utils.run_command(f"samtools view -bs {subsample_seed}{subsample_str} {input_bam} > {output_bam}")
                utils.run_command(f"samtools index -@ {threads} {output_bam}")
                pheno_bam_file_dct[pheno].append(new_bam)

    # Organize BAMs for each phenotype
    bam_list_files = []
    for p in subsample_pheno_dct:
        list_file = os.path.join(temp_dir, f"bam.{subsample_pheno_dct[p]}.list")
        with open(list_file, 'w') as wf:
            for bam in pheno_bam_file_dct[p]:
                wf.write(os.path.join(temp_dir, bam) + '\n')
            for sample in phenotype_dct[p]:
                if sample not in subsample_pheno_sample_dct[p]:
                    wf.write(os.path.join(temp_dir, f"{sample}.sorted.groupAdded.bam") + '\n')
        bam_list_files.append(list_file)
    
    logging.info("Finished BAM subsampling.")
    
    return bam_list_files


def call_variants(bam_list, output_dir, ref_seq, ploidy, threads):
    vcf_gz = os.path.join(output_dir, "combined.g.vcf.gz")
    vcf = os.path.join(output_dir, "combined.g.vcf")
    if not utils.file_exists(vcf):
        logging.info(f"Start variant calling.")
        utils.run_command(f"gatk HaplotypeCaller -R {ref_seq} -I {bam_list} "
                    f"-O {vcf_gz} --sample-ploidy {ploidy} --native-pair-hmm-threads {threads}")
        utils.run_command(f"gunzip -f {vcf_gz}")
        logging.info(f"Variant calling completed.")
    else:
        logging.info(f"VCF file exists. Skip variant calling.")
    return vcf


def annotate_variants(vcf, output_dir, snpeff_config):
    annotated_vcf = os.path.join(output_dir, "combined.g.annotated.vcf")
    if utils.file_exists(annotated_vcf):
        logging.info(f"Annotated VCF file already exists.")
        return annotated_vcf, True
    elif snpeff_config and utils.file_exists(vcf):
        logging.info(f"Start variant annotation.")
        utils.run_command(f"snpEff -config {snpeff_config} reference {vcf} > {annotated_vcf}")
        logging.info(f"Variant annotation completed.")
        return annotated_vcf, True
    else:
        logging.info(f"Skip variant annotation.")
        return vcf, False


def read_vcf(vcf_file, depth_cutoff, annotation_flag, filter_level):
    """Load variants from VCF file and filter variants by protein-coding impacts."""
    genotype_dct = {}
    with open(vcf_file, 'r') as rf:
        for line in rf:
            if line[:2]=='##':
                continue
            elif line[0]=='#': # Save sample names
                samples = line[:-1].split('\t')[9:]
            else:
                # Process each variant
                contig, pos, _, ref, alt = line.split('\t')[:5]
                alt_lst = alt.split(',')
                allele_lst = [ref]+alt_lst
                genotypes = line.split('\t')[9:]
                
                # Determine variant impact if filtering is enabled
                valid_genotype_lst = [0]
                best_gene_id = ""
                if annotation_flag:
                    level_dct = {l:i for i, l in enumerate(['MODIFIER', 'LOW', 'MODERATE', 'HIGH'])}
                    alt_dct = {al:i+1 for i, al in enumerate(alt_lst)}
                    ann_dct = {al:-1 for al in alt_lst}
                    gene_id_dct = {al: "" for al in allele_lst}
                    ann_field = next((f for f in line.split('\t')[7].split(';') if f.startswith('ANN=')), None)
                    if ann_field:
                        for ann in ann_field[4:].split(','):
                            cur_allele, cur_impact, cur_level, cur_gene_name, cur_gene_id = ann.split('|')[:5]
                            if level_dct[cur_level] > ann_dct[cur_allele]:
                                ann_dct[cur_allele] = level_dct[cur_level]
                                gene_id_dct[cur_allele] = f"{cur_gene_id}|{cur_level}|{cur_impact}"
                    
                    max_level = -1
                    for k in ann_dct:
                        level = ann_dct[k]
                        if level >= filter_level: # adjustable
                            valid_genotype_lst.append(alt_dct[k])
                        if level > max_level:
                            max_level = level
                            best_gene_id = gene_id_dct[k]
                    # Skip this variant if no allele meets the filter_level
                    if max_level < filter_level:
                        continue

                # Filter variants by depth and protein-coding impact
                variant_id = f"{contig}:{pos}|{best_gene_id}"
                genotype_dct[variant_id] = {}
                for k, g in enumerate(genotypes):
                    g_lst = []
                    depth_vals = g.split(':')[1].split(',')
                    for i, d in enumerate(depth_vals):
                        if d == '.':
                            continue
                        if int(d) >= depth_cutoff:
                            if not annotation_flag or filter_level == -1 or (i in valid_genotype_lst):
                                g_lst.append(i)
                    
                    if len(g_lst) == 0:
                        genotype_dct[variant_id][samples[k]] = -1
                    else:
                        genotype_dct[variant_id][samples[k]] = sum(g_lst)
    
    return genotype_dct


def run_pipeline(input_file, genome, genbank, output, force, mode, filter_level,
                 paired, temp, keep_temp, ploidy, min_depth, min_maf, min_sample_size, qcutoff, threads):
    logging.info(f"Pipeline started with {threads} threads")
    logging.info(f"Output will be stored at: {output}")
    os.makedirs(temp, exist_ok=True)
    logging.info(f"Temporary folder ready: {temp}")
    
    # 1. Reference preparation
    ref_seq, genes_lst, bt2_index, snpeff_config = prepare_reference(genome, genbank, temp, threads)
    
    # 2. Input and phenotype loading
    R1_list, R2_list, phenotype_dct = parse_input(input_file, paired)
    
    # 3. Mapping
    map_reads(R1_list, R2_list, paired, bt2_index, temp, threads)
    
    # 4. Compute coverate
    process_coverage(temp, threads)
    
    # 5. Analyze sequencing depth to determine if subsampling is needed
    subsample_stats_file = os.path.join(temp, "subsampling.stats")
    subsample_pheno_lst = analyze_sequencing_depths(temp, phenotype_dct, subsample_stats_file)
    
    # 6. Subsampling
    if subsample_pheno_lst:
        bam_list_files = subsample_bam_files(
            temp_dir=temp,
            output_dir=output,
            phenotype_dct=phenotype_dct,
            subsample_pheno_lst=subsample_pheno_lst,
            subsample_stats_file=subsample_stats_file,
            threads=threads
        )
    else:
        bam_list_files = []
        logging.info("No phenotypes require subsampling.")

    # Save a list of all original BAM files
    # if len(subsample_pheno_lst) < len(phenotype_dct):
    all_bam_file = os.path.join(temp, 'bam.all.list')
    with open(all_bam_file, 'w') as wf:
        for f in R1_list:
            sample = f.split('/')[-1]
            bam_file = f"{sample}.sorted.groupAdded.bam"
            wf.write(os.path.join(temp, bam_file) + '\n')
    bam_list_files.append(all_bam_file)

    # Run variant calling and GWAS for each BAM list file
    for bam_list_file in bam_list_files:
        subset_name = bam_list_file.split('.')[-2]
        cur_output = os.path.join(output, subset_name)
        os.makedirs(cur_output, exist_ok=True)

        logging.info(f"Processing subset: {subset_name}...")

        # 7. Variant calling
        vcf = call_variants(bam_list_file, cur_output, ref_seq, ploidy, threads)
        vcf, annotation_flag = annotate_variants(vcf, cur_output, snpeff_config)
        filter_level = int(filter_level)
        genotype_by_variant = read_vcf(vcf, min_depth, annotation_flag, filter_level)

        # 8. GWAS
        gwas.run_GWAS(genotype_by_variant, phenotype_dct, mode, min_maf, min_sample_size, qcutoff, cur_output)

    # Remove intermediate results
    if not keep_temp:
        shutil.rmtree(temp)
        logging.info("Temporary folder removed.")
    else:
        logging.info("Temporary folder kept for inspection.")