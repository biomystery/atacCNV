# List of scripts #

## Useage

```Shell
export PATH=$PATH:$(pwd)
```

## List of files 

* [bedToBdgBW.sh](./bedToBdgBw.sh): use `macs pileup` to generate pileup track (`bw`,`bdg`) from input reads (`bed`)  with reads extend half of window size in both direction. 
* [tagToFrag.sh](./tagToFrag.sh): covert `tagAlign` to fragments in `bed`, `sort.bdg` and `bw`  format. 
* [readsMasker.sh](./readsMasker.sh): filter reads (in `bed` format) by given mask regions. 
* [smoothBw.sh](./smoothBw.sh): smooth bigwig by sliding window (size - w and step - s)
* [callCNV.R](./callCNV.R): call CNV from treat adjusted by mapb 


## Typical Steps:

1. Get output from pipleine: `tagAlign`, `peaks` (openning regions).
2. (optional) Run `tagToFrag.sh`: get fragments `bed`, `bdg` and `bw`
3. Run `readsMasker.sh` to filter out reads or frags in peak regions and get background frags. 
4. Run `smoothBw.sh` to get average coverage in each window for both control or treatment. It will reorder the bed file, so that the it cannot preserve the pair info if input is in  `tagAlign` format. Please covert to fragment bed first. 
5. Run `callCNV.sh` to get CNV regions, segment plots 

## Call Wrapper: 
* Input: reads (`tagAlign`), peak regions (`bed` or `narrowPeaks`), mapb track (`bw`), sliding window (`stepsize`,`windowsize`)
* Output: CNV regions (`bed`), log2CovRaito plots with segments (`pngs` and `tab`), background reads track (`bw`) 


