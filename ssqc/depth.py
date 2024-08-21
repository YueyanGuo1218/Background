import pysam
from typing import NamedTuple
import numpy as np

def spikiness(chr_list, ref_list, ref_lengths, watson_starts, crick_starts):
    result = np.zeros(len(chr_list))
    count = 1
    for i in calculate_depth(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, [], 100_000, 75):
        result[count-1] = np.sum(np.abs(np.diff(i))) / np.sum(i)
        count += 1
    return np.mean(result)

def half_depth_proportion(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, total_depth, mask_dir = None, window_size = 1_000_000):
    def sum_continuous_ones(arr):
        # shift arr right
        r = np.roll(arr, 1)
        r[0] = 0

        # shift arr rleft
        l = np.roll(arr, -1)
        l[-1] = 0

        starts = np.where(arr - r == 1)[0] # if an element is greater than the element on its left, it is a start
        ends = np.where(arr - l == 1)[0]

        # Replace sequences of 1s with the count of 1s
        for start, end in zip(starts, ends):
            count = end - start + 1
            arr[start:end+1] = 0  # Reset the range to 0
            arr[start] = count  # Set the start of the range to the count

        return arr
    total_reads = total_depth
    total_length = 0
    for chr in chr_list:
        if chr not in ref_list:
            continue
        total_length += ref_lengths[chr]
    proportion_of_single_window = window_size / total_length
    

    result = np.zeros(len(chr_list))
    count = 1
    for i in calculate_depth(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, [], window_size, 75):
        i = i / (total_reads * proportion_of_single_window) # normalized depth
        half_depth_array = np.where(np.abs(i-0.5) < 0.1, 1, 0)
        #half_depth_ratio = np.sum(half_depth_array) / len(half_depth_array)
        #print(half_depth_array)
        half_depth_ratio = np.square(sum_continuous_ones(half_depth_array))
        #print(half_depth_array)
        half_depth_ratio = np.sum(half_depth_array) / len(half_depth_array)
        result[count-1] = half_depth_ratio
        #print(f"chr{count}: {half_depth_ratio}")
        count += 1
        #print(result)
    return np.nanmean(result)

def etp(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, mask_dir = None, window_size = 1_000_000):
    def entropy_formula(arr: np.ndarray):
        arr = arr / np.sum(arr)
        return -np.nansum(arr * np.log2(arr))

    result = np.zeros(len(chr_list))
    count = 1

    for i in calculate_depth(chr_list, ref_list, ref_lengths, watson_starts, crick_starts, [], 100_000, 75):
        result[count-1] = entropy_formula(i)
        count += 1
    return np.nanmean(result)

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

def calculate_depth(chr_list,
                    ref_name, ref_length,
                    watson_starts, crick_starts,
                    masked_regions, bin_size, read_length):
    for chr in chr_list:
        if chr not in ref_name:
            print(f"{chr} not in reference (bam header section)")
            continue
        chr_id = ref_name.index(chr)
        chr_length = ref_length[chr]
        yield calculate_depth_chr(chr_id, chr_length, watson_starts, crick_starts, masked_regions, bin_size, read_length)

def calculate_depth_chr(chr_id, chr_length, watson_starts, crick_starts, masked_regions, bin_size, read_length):
    size = array_size(chr_length, bin_size)
    w_count = np.zeros(size)
    c_count = np.zeros(size)

    for read in watson_starts:
        if read.chr_id == chr_id:
            w_count[read.start] += 1
        else:
            break #assuming the bam file is sorted
            #continue
    w_count = sum_by_bin(w_count, bin_size)

    for read in crick_starts:
        if read.chr_id == chr_id:
            c_count[read.start] += 1
        else:
            #break #assuming the bam file is sorted
            continue
    c_count = sum_by_bin(c_count, bin_size)

    depth = w_count + c_count

    #print(f"working on chr{chr_id+1}")
    #print(depth)
    #print(relative_depth)
    #print(relative_depth > 0.6)

    return depth
