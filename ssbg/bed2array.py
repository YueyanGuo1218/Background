import pysam
from typing import NamedTuple
import numpy as np
from typing import NamedTuple, List
from pathlib import Path

Region = NamedTuple('Region', [('chr', str), ('start', int), ('end', int)])

def read_bedfile(bed_file):
    ''' return a list of regions '''
    regions = []
    with open(bed_file, 'r') as f:
        for line in f:
            chr, start, end = line.strip().split()
            start = int(start)
            end = int(end)
            regions.append(Region(chr, start, end))
    return regions

def get_reference(bam_file):
    ''' return reference sequence name and length from a bamfile'''

    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get reference seqs name and length
        ref_list = bam.references
        ref_lengths = {}
        for reference in bam.header['SQ']:
            chrom_name = reference['SN']
            chr_length = reference['LN']
            ref_lengths[chrom_name] = chr_length
    return ref_list, ref_lengths

def array_size(length, bin_size):
    bin_count = length // bin_size
    if length % bin_size != 0:
        bin_count += 1
    return bin_count * bin_size

def sum_by_bin(array, bin_size):
    size = array_size(len(array), bin_size)
    reshaped_array = array.reshape((size // bin_size, bin_size))
    reshaped_array = np.sum(reshaped_array, axis=1)
    return reshaped_array

def bed2array(mask_regions: List[Region], reference_list, reference_length, window_size = 1_000_000) -> List[np.array]:
    ''' return a list of arrays with maskd window labeld 0 '''

    chromosome_list = [f"chr{i}" for i in range(1, 23)]

    for chr in chromosome_list:
        if chr not in reference_list:
            continue

        chr_length = reference_length[chr]
        size = array_size(chr_length, window_size)
        mask_array = np.ones(size)

        for region in mask_regions:
            if region.chr == chr:
                mask_array[region.start:region.end] = 0

        mask_array = sum_by_bin(mask_array, window_size)
        mask_array = np.where(mask_array < window_size, 0, 1)
        yield mask_array

def generate_array(bedfile: str, bamfile: str, chr_list: List[str], window_size: int, output_dir: str):
    ''' save the array file (.npy) for every chromosome in chr_list into the specified output directory'''
    regions = read_bedfile(bedfile)
    ref_list, ref_length = get_reference(bamfile)

    results = list(bed2array(regions, ref_list, ref_length, window_size)) # a list of arrays with masked window labeled 0

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    for chr in chr_list:
        mask = results[chr_list.index(chr)]
        np.save(output_path / f"{chr}.npy", mask)
