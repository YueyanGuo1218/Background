import pysam
import numpy as np
from pathlib import Path
from .background import process_bam_file, calculate_background
#from .bed2array import bed2array

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Process BAM files and BED masks.')
    parser.add_argument('input_bam', type=str, help='Input BAM file')
    parser.add_argument('mask_bed', type=str, help='Mask BED file')
    args = parser.parse_args()
    
    # Your code to process the BAM and BED files
    print(f'Processing {args.input_bam} with mask {args.mask_bed}')
    if args.input_bam == "-":
        bam_file = None
    else:
        bam_file = args.input_bam

    ref_list, ref_lengths, watson_starts, crick_starts = process_bam_file(bam_file)

    mask_dir = Path(args.mask_bed).parent
    cache_folder_name = "ssbg_cache"
    cache_path = mask_dir / cache_folder_name
    #look for cache in the dir of bed file
    if cache_path.is_dir():
        print(f"Using cached files in {cache_path}")
    else:
        print(f"Generating cache for {args.mask_bed}")

    chr_list = [f"chr{i}" for i in range(1,23)]
    background = calculate_background(chr_list, 
                         ref_list, ref_lengths, watson_starts, crick_starts, 
                         masked_regions = None, bin_size = 1_000_000, read_length = 75)
    print(background)


if __name__ == '__main__':
    main()
