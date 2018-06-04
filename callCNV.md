# callCNV.R

## Usage

```Bash
$ Rscript callCNV.R -h

Usage: callCNV.R [options]


Options:
        -c CHARACTER, --ctrl=CHARACTER
                ctrl bed file [default= ./data/hg38_75bp/GRCh38_75.avg_cov.w1000K_s500K.bed]

        -t CHARACTER, --treat=CHARACTER
                treat bed file

        -g CHARACTER, --genome=CHARACTER
                genome size file [default= ../data/hg38.chrom.sizes]

        -s CHARACTER, --sample=CHARACTER
                sample name

        -u DOUBLE, --upper=DOUBLE
                upper boundary in trimming log2Ratio [default= 2]

        -m DOUBLE, --mapb_th=DOUBLE
                mapbility threshold for regions [default= 0.75]

        -o CHARACTER, --output_dir=CHARACTER
                output dir [default= ./data/output/]

        -l DOUBLE, --lower=DOUBLE
                lower boundary in trimming log2Ratio [default= -1.5]

        -h, --help
                Show this help message and exit

```
