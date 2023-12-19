configfile:  '../config/config.yaml'  #'/home/people/lasdan/lasdan/benchmark_workflow/config/config.yaml'
SCRIPTS = '/home/people/lasdan/lasdan/vamb_versions/vamb/src/'
CUSTOM_SCRIPTS =  '../workflow/scripts/' #'/home/people/lasdan/lasdan/benchmark_workflow/workflow/scripts/'    

VAMBCONDAENV = 'vamb' 
SAMTOOLSENV = "samtools/1.17"
MINIMAP2ENV = "minimap2/2.6"

import re
import os
import sys
import numpy as np
SNAKEDIR = os.path.dirname(workflow.snakefile)

sys.path.append(os.path.join(SNAKEDIR, 'scripts'))

def get_config(name, default, regex):
    res = config.get(name, default).strip()
    m = re.match(regex, res) 
    if m is None:
        raise ValueError(
            f"Config option \"{name}\" is \"{res}\", but must conform to regex \"{regex}\"")
    return res

# set configurations
CONFIG_PATH = config['config_path']
SAMPLE = config['sample']
CONTIGS = get_config("contigs", CONFIG_PATH + SAMPLE + "/contigs.txt", r".*") # each line is a contigs path from a given sample
TIMES_RUN_VAMB = config['times_run_vamb']

SAMPLE_DATA =  CONFIG_PATH + SAMPLE + "/samples2data.tsv"
# TODO remove 
config['contigs'] = CONTIGS
# print('SAMPLE_DATA:', SAMPLE_DATA)

INDEX_SIZE = config['index_size'] 
MIN_BIN_SIZE = config['min_bin_size']
MIN_IDENTITY = config['min_identity'] 
MM_MEM = config['minimap']['minimap_mem']
MM_PPN = config['minimap']['minimap_ppn']
AVAMB_MEM = config['vamb']['avamb_mem']
AVAMB_PPN = config['vamb']['avamb_ppn']


MIN_COMP = config['min_comp']
MAX_CONT = config['max_cont']



# parse if GPUs is needed #
AVAMB_PPN = AVAMB_PPN 
AVAMB_GPUS = config['vamb']['avamb_gpus']
CUDA = AVAMB_GPUS > 0

## read in sample information ##

IDS = []
sample2path = {}
fh_in = open(SAMPLE_DATA, 'r')
LONG_READS = None
one_samples = False
two_samples = False
head_path = False
path = ''
for line in fh_in:
    line = line.rstrip()
    fields = line.split()

    ### Basic Reading of file
    # Checking for mistakes in the files
    if two_samples and len(fields) == 2:
        raise Exception('both 1 and 2 samples')
    if one_samples and len(fields) == 3:
        raise Exception('both 1 and 2 samples')
#    if head_path and len(fields) == 1:
#        raise Exception('1 sample line after header')  

    # Adding the data for entries with 1
    if len(fields) == 3:
        sample2path[fields[0]] = [path + fields[1], path + fields[2]]
        two_samples = True
        # print('Assuming paired reads, since given 2 fastq sequences per sample:')
        LONG_READS = False
        # print('LONG_READS=',LONG_READS)

        IDS.append(fields[0])

    
    elif len(fields) == 2:
        sample2path[fields[0]] = [path + fields[1]]
        one_samples = True
        # print('Assuming long-reads, since given 1 fastq sequence per sample')
        LONG_READS = True
        # print('LONG_READS=',LONG_READS)

        IDS.append(fields[0])

    elif len(fields) == 1:
        path = fields[0]
        head_path = True



#
# read in list of per-sample assemblies
CONTIGS = CONFIG_PATH + SAMPLE + "/contigs.txt"
contigs_list = []
fh_in = open(CONTIGS, 'r')
for line in fh_in:
    line = line.rstrip()
    contigs_list.append(line)



############
vamb_cutoff_list = [500, 1000, 1500, 2000]

rule all:
    input:
        outdir_avamb=expand('cutoff_{vamb_cutoff}/{vamb_runs}_vamb/vae_clusters_split.tsv',vamb_runs = list(range(1, TIMES_RUN_VAMB + 1)), vamb_cutoff = vamb_cutoff_list)

# Filter contigs for 2000bp and rename them to conform with the multi-split workflow 
rule cat_contigs:
    input:
        contigs_list
    output:
        "contigs.flt.fna.gz"
    params:
        path=SCRIPTS + 'concatenate.py',
        walltime="864000",
        nodes="1",
        ppn="1",
    resources:
        mem="5GB"
    threads:
        1
    conda: 
        VAMBCONDAENV
    # log:
        # o = os.path.join(OUTDIR,"log/contigs/catcontigs.o"),
        # e = os.path.join(OUTDIR,"log/contigs/catcontigs.e")
    shell: 
        "python {params.path} {output} {input} -m 0" 

# Index resulting contig-file with minimap2
rule index:
    input:
        contigs = "contigs.flt.fna.gz"
    output:
        mmi = "contigs.flt.mmi"
    params:
        walltime="864000",
        nodes="1",
        ppn="1"
    resources:
        mem="90GB"
    threads:
        1
    group:
        '2_index_and_run_minimap2'
    log:
        out_ind ="log/contigs/index.log",
        o = "log/contigs/index.o",
        e = "log/contigs/index.e",
 
    envmodules:
        'tools',
        config['moduleenvs']['minimap2']
    #conda: 
    #    "envs/minimap2.yaml"
    shell:
        "minimap2 -I {INDEX_SIZE} -d {output} {input} 2> {log.out_ind}"

# This rule creates a SAM header from a FASTA file.
# We need it because minimap2 for truly unknowable reasons will write
# SAM headers INTERSPERSED in the output SAM file, making it unparseable.
# To work around this mind-boggling bug, we remove all header lines from
# minimap2's SAM output by grepping, then re-add the header created in this
# rule.
rule dict:
    input:
        contigs = "contigs.flt.fna.gz",
    output:
        dict = "contigs.flt.dict"
    params:
        walltime="864000",
        nodes="1",
        ppn="1"
    resources:
        mem="10GB"
    threads:
        1
    log:
        out_dict= "log/contigs/dict.log",
        o = "log/contigs/dict.o",
        e = "log/contigs/dict.e"
    envmodules:
        'tools',
        config['moduleenvs']['samtools']
    group:
        '2_index_and_run_minimap2'

    shell:
        "samtools dict {input} | cut -f1-3 > {output} 2> {log.out_dict}"

# Generate bam files 
rule minimap:
    input:
        fq = lambda wildcards: sample2path[wildcards.sample],
        mmi ="contigs.flt.mmi",
        dict = "contigs.flt.dict"
    output:
        bam = temp("mapped/{sample}.bam")
    params:
        walltime="864000",
        nodes="1",
        ppn=MM_PPN,
        long_or_short_read = 'map-pb -L' if LONG_READS else 'sr',
    resources:
        mem=MM_MEM
    threads:
        int(MM_PPN)
    log:
        out_minimap = "log/map/{sample}.minimap.log",
        o = "log/map/{sample}.minimap.o",
        e = "log/map/{sample}.minimap.e",
    envmodules:
        'tools',
        config['moduleenvs']['minimap2'],
        config['moduleenvs']['samtools']
    group:
        '2_index_and_run_minimap2'
    shell:
        # See comment over rule "dict" to understand what happens here
        "minimap2 -t {threads} -ax {params.long_or_short_read} {input.mmi} {input.fq} -N 5"
        " | grep -v '^@'"
        " | cat {input.dict} - "
        " | samtools view -F 3584 -b - " # supplementary, duplicate read, fail QC check
        " > {output.bam} 2> {log.out_minimap}"

# Sort bam files
rule sort:
    input:
        "mapped/{sample}.bam",
    output:
        "mapped/{sample}.sort.bam",
    params:
        walltime="864000",
        nodes="1",
        ppn="2",
        prefix="mapped/tmp.{sample}"
    resources:
        mem="15GB"
    threads:
        2
    log:
        out_sort = "log/map/{sample}.sort.log",
        o = "log/map/{sample}.sort.o",
        e = "log/map/{sample}.sort.e",
    envmodules:
       'tools',
        config['moduleenvs']['samtools']
        #conda:
    #    '/mnt/c/Users/Admin/Desktop/benchmark_workflow/workflow/envs/samtools.yaml'
    shell:
        "samtools sort {input} -T {params.prefix} --threads 1 -m 3G -o {output} 2> {log.out_sort}"

rule run_avamb:
    input:
        contigs="contigs.flt.fna.gz",
        bam_files = expand("mapped/{sample}.sort.bam", sample = IDS)
    output:
        outdir_avamb=directory('cutoff_{vamb_cutoff}/{vamb_runs}_vamb'),
        file = 'cutoff_{vamb_cutoff}/{vamb_runs}_vamb/vae_clusters_split.tsv'
    params:
        walltime="86400",
        nodes="1",
    resources:
        mem=AVAMB_MEM
    threads:
        int(8)
    conda:
        VAMBCONDAENV
    log:
        vamb_out="tmp/{vamb_cutoff}/{vamb_runs}_vamb_finished.log",
        o=os.path.join('log','cutoff_{vamb_cutoff}/{vamb_runs}_vamb.out'),
        e=os.path.join('log','cutoff_{vamb_cutoff}/{vamb_runs}_vamb.err')
    shell:
        "rm -rf {output.outdir_avamb} || true;"
        "vamb bin default -z 0.95 --outdir {output.outdir_avamb} --fasta {input.contigs} "
        "--bamfiles {input.bam_files} -m {wildcards.vamb_cutoff} "



