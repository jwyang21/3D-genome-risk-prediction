from pathlib import Path
from os import path

REFERENCE_DIR = Path(config['reference']['dir'])
REFERENCE_NAME = config['reference']['name']
INDEX_DIR = Path(config['bwa_index_dir'])

def bwa_mem_input(wildcards):

    ret = {}
    ret['reference'] = str(INDEX_DIR / f'{REFERENCE_NAME}.bwt')
    ret['reads'] = [
        str(data_dir / f'{wildcards.run}.read1.fastq.gz'),
        str(data_dir / f'{wildcards.run}.read2.fastq.gz'),
    ]
    return ret

c = config['bwa_index']
rule bwa_index:
    input:
        # Required input. Reference genome fasta file.
        config['reference']['fasta'],
    output:
        # Required output. BWA-indexed reference genome files.
        INDEX_DIR / (REFERENCE_NAME + '.amb'),
        INDEX_DIR / (REFERENCE_NAME + '.ann'),
        INDEX_DIR / (REFERENCE_NAME + '.bwt'),
        INDEX_DIR / (REFERENCE_NAME + '.pac'),
        INDEX_DIR / (REFERENCE_NAME + '.sa'),
    params:
        extra = c['extra'],
        # Note that the default algorithm for this wrapper is 'bwtsw'.
        # The other option is 'is', but please be warned that this algorithm doesn't work
        # with genomes longer than 2GB.
        # Default: 'bwtsw',
        a = c['a'],
        # Block size for the bwtsw algorithm (effective with -a bwtsw).
        # Default: False [10000000]
        b = c['b'],
        # Index files named as <in.fasta>.64.* instead of <in.fasta>.*
        _6 = c['_6'],
    threads: config['threads']['bwa_index']
    log: f'logs/bwa_index/{REFERENCE_NAME}.log'
    benchmark: f'benchmarks/bwa_index/{REFERENCE_NAME}.log'
    wrapper:
        'http://dohlee-bio.info:9193/bwa/index'

c = config['bwa_mem']
rule bwa_mem:
    input: unpack(bwa_mem_input)
    output:
        # BAM output or SAM output is both allowed.
        # Note that BAM output will be automatically detected by its file extension,
        # and SAM output (which is bwa mem default) will be piped through `samtools view`
        # to convert SAM to BAM.
        result_dir / '01_bwa' / '{run}.bam'
    params:
        extra = c['extra'],
        # Minimum seed length.
        # Default: 19
        k = c['k'],
        # Band width for banded alignment.
        # Default: 100
        w = c['w'],
        # Off-diagonal X-dropoff.
        # Default: 100
        d = c['d'],
        # Look for internal seeds inside a seed longer than {-k} * FLOAT
        # Default: 1.5
        r = c['r'],
        # Seed occurrence for the 3rd round seeding.
        # Default: 20
        y = c['y'],
        # Skip seeds with more than INT occurrences.
        # Default: 500
        c = c['c'],
        # Drop chains shorter than FLOAT fraction of the logest overlapping chain.
        # Default: 0.5
        D = c['D'],
        # Discard a chain if seeded bases shorter than INT.
        # Default: 0
        W = c['W'],
        # Perform at most INT rounds of mate rescues for each read.
        # Default: 50
        m = c['m'],
        # Skip mate rescue.
        # Default: False
        S = c['S'],
        # Skip pairing; mate rescue performed unless -S also in use
        # Default: False
        P = c['P'],
        # Score for a sequence match, which scales options -TdBOELU unless overridden.
        # Default: 1
        A = c['A'],
        # Penalty for a mismatch.
        # Default: 4
        B = c['B'],
        # Gap open penalties for deletions and insertions.
        # Default: 6,6
        O = c['O'],
        # Gap extension penalty; a gap of size k cost '{-O} + {-E}*k'.
        # Default: 1,1
        E = c['E'],
        # Penalty for 5'- and 3'-end clipping.
        # Default: 5,5
        L = c['L'],
        # Penalty for an unpaired read pair.
        # Default: 17
        U = c['U'],
        # Read type. Setting -x changes multiple parameters unless overridden [null]
        # pacbio: -k17 -W40 -r10 -A1 -B1 -O1 -E1 -L0 (PacBio reads to ref)
        # ont2d: -k14 -W20 -r10 -A1 -B1 -O1 -E1 -L0 (Oxford Nanopore 2D-reads to ref)
        # intractg: -B9 -O16 -L5 (intra-species contigs to ref)
        # Default: False
        x = c['x'],
        # Read group header line such as '@RG\tID:foo\tSM:bar'
        # Default: False
        # NOTE: You should check the platform information of the read data!
        R = c['R'],
        # Insert STR to header if it starts with @; or insert lines in FILE.
        # Default: False
        H = c['H'],
        # Treat ALT contigs as part of the primary assembly. (i.e. ignore <idxbase>.alt file)
        # Default: False
        j = c['j'],
        # For split alignment, take the alignment with the smallest coordinate as primary.
        # Default: False
        _5 = c['_5'],
        # Dont't modify mapQ of supplementary alignments.
        # Default: False
        q = c['q'],
        # Process INT input bases in each batch regardless of nThreads (for reproducibility).
        # Default: False.
        K = c['K'],
        # Verbosity level: 1=error, 2=warning, 3=message, 4+=debugging
        # Default: 3
        v = c['v'],
        # Minimum score to output.
        # Default: 30
        T = c['T'],
        # If there are <INT hits with score > 80% of the max score, output all in XA.
        # Default: 5,200
        h = c['h'],
        # Output all alignments for SE or unpaired PE.
        # Default: False
        a = c['a'],
        # Append FASTA/FASTQ comment to SAM output.
        # Default: False
        C = c['C'],
        # Output the reference FASTA header in the XR tag.
        # Default: False
        V = c['V'],
        # Use soft clipping for supplementary alignments.
        # Default: False
        Y = c['Y'],
        # Mark shorter split hits as secondary.
        # NOTE: You may need this if you use GATK downstream.
        # Default: False
        M = c['M'],
        # Specify the mean, standard deviation (10% of the mean if absent), max
        # (4 sigma from the mean if absent) and min of the insert size distribution.
        # FR orientation only.
        # Default: False (inferred)
        I = c['I'],
        index = lambda wc, input: os.path.splitext(input.reference)[0],
    threads: config['threads']['bwa']
    log: 'logs/bwa_mem/{run}.log'
    benchmark: 'benchmarks/bwa_mem/{run}.benchmark'
    shell:
        'bwa mem -t {threads} -SP5M {params.index} {input.reads} | samtools view -Shb - > {output}'
