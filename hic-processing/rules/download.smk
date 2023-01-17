ASCP_BIN = config['ascp_bin']
ASCP_KEY = config['ascp_key']

rule prefetch_accession:
    output:
        temp('{run}.sra')
    resources:
        network = 1
    shell:
        'prefetch --ascp-path "{}|{}" -v {{wildcards.run}} --max-size 1000000000 && mv {{wildcards.run}}/{{wildcards.run}}.sra . && rm -r {{wildcards.run}}'.format(ASCP_BIN, ASCP_KEY)

rule parallel_fastq_dump_paired:
    input:
        # Required input. Recommend using wildcards for sample names,
        # e.g. {sample,SRR[0-9]+}
        '{sample}.sra'
    output:
        # Required output.
        data_dir / '{sample}.read1.fastq.gz',
        data_dir / '{sample}.read2.fastq.gz',
        # temp('{sample}_pass.fastq.gz')
    params:
        # Optional parameters. Omit if unused.
        extra = '--tmpdir .'
    threads: 4
    wrapper:
        'http://dohlee-bio.info:9193/parallel-fastq-dump'

