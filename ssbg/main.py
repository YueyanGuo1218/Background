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
    parser.add_argument('--output', type=str, default = ".", help='Output directory')
    parser.add_argument('--window_size', type=int, default = 1_000_000, help='Window size used in calculation')
    args = parser.parse_args()

    chr_list = [f"chr{i}" for i in range(1,23)]

    print(f'Processing {args.input_bam}')

    # Process input_bam arguement
    if args.input_bam == "-":
        bam_file = None     # calculate_background will use stdin
    else:
        bam_file = args.input_bam

    ref_list, ref_lengths, watson_starts, crick_starts = process_bam_file(bam_file)

    # Process output arguement
    if args.output == ".":
        output_path = Path.cwd()
    else:
        output_path = Path(args.output)
        output_path.mkdir(parents=True, exist_ok=True)

    # Process mask argument
    if args.mask is None:
        cache_path = None
    else:
        bed_file = Path(args.mask)
        mask_dir = bed_file.parent
        cache_folder_name = f"ssbg_cache_{bed_file.stem}"
        cache_path = mask_dir / cache_folder_name

        #look for cache in the dir of bed file
        if cache_path.is_dir():
            print(f"Using cached mask arrays in {cache_path}")
        else:
            generate_array(bed_file, bam_file, chr_list, args.window_size, cache_path)
            print(f"Generating mask arrays for {args.mask}")


    background = calculate_background(chr_list,
                         ref_list, ref_lengths, watson_starts, crick_starts,
                         mask_dir = cache_path, bin_size = 1_000_000, read_length = 75)
    print(f"Background over whole genome: {background}")

    if bam_file is None:
        with open(output_path / f"stdin_background.txt", "w") as f:
            f.write(f"{background}")
    else:
        with open(output_path / f"{Path(args.input_bam).stem}_background.txt", "w") as f:
            f.write(f"{background}")


if __name__ == '__main__':
    main()
