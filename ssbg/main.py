import pysam
import numpy as np
from pathlib import Path
from .background import process_bam_file, calculate_background
from .bed2array import generate_array

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Process BAM files and BED masks.')
    parser.add_argument('input_bam', type=str, help='Input BAM file')
    parser.add_argument('--mask', type=str, default = None, help='Mask BED file')
    parser.add_argument('--window_size', type=int, default = 1_000_000, help='Window Size used in calculation')
    args = parser.parse_args()

    chr_list = [f"chr{i}" for i in range(1,23)]

    print(f'Processing {args.input_bam}')

    if args.input_bam == "-":
        bam_file = None     # calculate_background will use stdin
    else:
        bam_file = args.input_bam

    ref_list, ref_lengths, watson_starts, crick_starts = process_bam_file(bam_file)

    # Process mask argument
    if args.mask is None:
        cache_path = None
    else:
        bed_file = Path(args.mask)
        mask_dir = bed_file.parent
        cache_folder_name = "ssbg_cache"
        cache_path = mask_dir / cache_folder_name

        #look for cache in the dir of bed file
        if cache_path.is_dir():
            print(f"Using cached files in {cache_path}")
        else:
            generate_array(bed_file, bam_file, chr_list, args.window_size, cache_path)
            print(f"Generating cache for {args.mask}")


    background = calculate_background(chr_list,
                         ref_list, ref_lengths, watson_starts, crick_starts,
                         mask_dir = cache_path, bin_size = 1_000_000, read_length = 75)
    print(background)


if __name__ == '__main__':
    main()
