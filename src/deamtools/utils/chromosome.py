import pysam


def get_chrom_sizes_from_bam(bam: pysam.Samfile) -> dict[str, int]:
    """
    Extract chromosome sizes from a BAM file.

    Parameters
    ----------
    bam : pysam.Samfile
        An open BAM file object from which chromosome names and lengths are read.

    Returns
    -------
    dict[str, int]
        A dictionary mapping chromosome names to their sizes (length - 1).
    """
    chromosome = list(bam.references)
    lengths = list(bam.lengths)
    chrom_sizes = dict(zip(chromosome, lengths))
    return chrom_sizes


def get_chrom_sizes_from_file(chrom_size_file: str) -> dict[str, int]:
    """
    Read chromosome sizes from a tab-delimited chromosome sizes file.

    Parameters
    ----------
    chrom_size_file : str
        Path to a tab-delimited file where each line contains a chromosome name
        and its size (e.g., UCSC-style chrom.sizes format).

    Returns
    -------
    dict[str, int]
        A dictionary mapping chromosome names to their sizes.
    """
    chrom_sizes = {}
    with open(chrom_size_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            size = int(parts[1])
            chrom_sizes[chrom] = size
    return chrom_sizes