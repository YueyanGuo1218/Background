import pysam
from typing import NamedTuple
import numpy as np
import sys

Read = NamedTuple('Read', [('chr_id', int), ('start', int)])

def process_bam_file(bam_file = None):
    watson_list = []
    crick_list = []
    if bam_file is None:
        bam_file = sys.stdin.buffer
    # Open BAM file from stdin (binary mode)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get reference seqs name and length
        ref_list = bam.references
        ref_lengths = {}
        for reference in bam.header['SQ']:
            chrom_name = reference['SN']
            chr_length = reference['LN']
            ref_lengths[chrom_name] = chr_length

        # Iterate over reads in BAM file
        for read in bam:
            flag = read.flag
            chr_id = read.reference_id
            start = read.reference_start
            read = Read(chr_id, start)
            # Check flag conditions
            if flag in [83, 163, 83+1024, 163+1024]:
                watson_list.append(read)
            elif flag in [99, 147, 99+1024, 147+1024]:
                crick_list.append(read)
            else:
                continue  # Ignore other flags

    return ref_list, ref_lengths, watson_list, crick_list

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

def strand_genotype(arr):
    transformed_arr = np.where(np.abs(arr) < 0.2, 0,
                    np.where(np.abs(arr) < 0.8, np.sign(arr)* 0.5,
                        np.sign(arr)))
    return transformed_arr

def calculate_background_chr(chr_id, chr_length, watson_starts, crick_starts, masked_regions, bin_size, read_length):
    size = array_size(chr_length, bin_size)
    w_count = np.zeros(array_size(chr_length, bin_size))
    c_count = np.zeros(array_size(chr_length, bin_size))

    for read in watson_starts:
        if read.chr_id == chr_id:
            w_count[read.start] += 1
        else:
            #break #assuming the bam file is sorted
            continue
    w_count = sum_by_bin(w_count, bin_size)
    for read in crick_starts:
        if read.chr_id == chr_id:
            c_count[read.start] += 1
        else:
            #break #assuming the bam file is sorted
            continue
    c_count = sum_by_bin(c_count, bin_size)
    strand = strand_genotype( (w_count - c_count) / (w_count + c_count) )

    print(f"working on chr{chr_id+1}")
    print(strand)

    for region in masked_regions:
        region_chr_id = ref_list.index(region.chr)
        if region_chr_id == chr_id:
            strand[region.start:region.end] = np.nan # label masked region as no reads

    ww_windows = (strand == 1)
    cc_windows = (strand == -1)

    ww_window_counts = np.count_nonzero(ww_windows)
    cc_window_counts = np.count_nonzero(cc_windows)

    if ww_window_counts == 0 and cc_window_counts == 0:
        return np.nan, 0
    elif ww_window_counts == 0:
        cc_background = np.average(w_count[cc_windows]/c_count[cc_windows])
        return cc_background, cc_window_counts
    elif cc_window_counts == 0:
        ww_background = np.average(c_count[ww_windows]/w_count[ww_windows])
        return ww_background, ww_window_counts


    print("ww & cc exist") # debug
    # if both cc ww regions exists, return weighted average

    ww_background = np.average(c_count[ww_windows]/w_count[ww_windows])
    cc_background = np.average(w_count[cc_windows]/c_count[cc_windows])

    ww_weight = ww_window_counts / (ww_window_counts + cc_window_counts)
    cc_weight = cc_window_counts / (ww_window_counts + cc_window_counts)

    return ww_weight * ww_background + cc_weight * cc_background, ww_window_counts + cc_window_counts

def calculate_background(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, masked_regions = None, bin_size = 1_000_000, read_length = 75):
    if masked_regions is None:
        masked_regions = []
    background = np.empty(len(chr_list))
    chr_weight = np.zeros_like(background)
    background.fill(np.nan)

    for chr in chr_list:
        if chr not in ref_list:
            print(f"Warning: {chr} not found in BAM file")
            continue
        chr_id = ref_list.index(chr)
        chr_length = ref_lengths[chr]
        background[chr_list.index(chr)], chr_weight[chr_list.index(chr)] = calculate_background_chr(chr_id,
                                                                                                    chr_length,
                                                                                                    watson_starts,
                                                                                                    crick_starts,
                                                                                                    masked_regions,
                                                                                                    bin_size,
                                                                                                    read_length)
        print(background[chr_list.index(chr)]) # debug

    #return background, chr_weight

    chr_weight = chr_weight / np.sum(chr_weight)
    bg_result = np.average(background[~np.isnan(background)], weights=chr_weight[~np.isnan(background)])
    print(bg_result)
    return bg_result
    
