import subprocess

#### SETTINGS ####
# TARGET_PARAMETER = [2000]
# PARAMETER_NAME = 'para' # VIRKER IKKE ENDU skal aendres de andre steder
# SAMPLES =  ['A'] #['Gastrointestinal', 'Skin', 'Urogenital'] # ['errorfree', 'sludge', 'human_longread']
SAMPLE_DATA_PATH = '../sample_data/'
IN_LOOP_SNAKEFILE = './vamb_cutoff_in_loop_snakefile.py'
TIMES_RUN_VAMB = 5 #! Change to 5 
# TARGET_PARAMETER = [500, 1000, 1500, 2000]
SAMPLES =  ['Airways', 'Gastrointestinal', 'Skin', 'Urogenital', 'Oral'] # ['errorfree', 'sludge', 'human_longread'

#### CODE ####
TIMES_RUN_VAMB_LIST = list(range(1, TIMES_RUN_VAMB + 1)) # list for expand in rule all input and run_inloop output

rule all:
    input:
         directory(expand('{sample}', sample = SAMPLES))
    
rule run_inloop:
    output: 
        directory('{sample}')
    params:
        config_path = SAMPLE_DATA_PATH,
        in_loop_snakefile = IN_LOOP_SNAKEFILE,
        times_run_vamb = TIMES_RUN_VAMB
    benchmark: 
        'benchmark/{sample}'
    shell:
        'mkdir  {wildcards.sample} -p;'
        'source activate snakemake;'
        'snakemake -c4 -d {wildcards.sample} --snakefile {params.in_loop_snakefile} ' 
        '--config config_path={params.config_path} '
        'sample={wildcards.sample} '
        'times_run_vamb={params.times_run_vamb} '
        '--use-envmodules --use-conda '
        # '-np '
        # '--rerun-incomplete'
        '2> output/{wildcards.sample}_snakemake.e '
        '> output/{wildcards.sample}_snakemake.o '



         


         

