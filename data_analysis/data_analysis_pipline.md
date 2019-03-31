# MRMix data analysis pipeline

### Preprocessing summary level data
Code for preprocessing and links to summary level data are provided [here](https://github.com/gqi/MRMix/blob/master/data_analysis/preprocess_summary_stats.R). 

To reproduce the results, first download the following files:
* A PLINK `bim` file for HapMap 3 SNPs with MAF>0.05 1000 Genomes Phase 3 European sample (around 1.07 million SNPs). Download [here](https://jh.box.com/s/qkp4h90nm6tp2qj82v5moxqs0no0gvwa).
* Linkage disequilibrium (LD) coefficients calculated from 1000G Phase 3 European sample, with ~1.2 million HapMap3 SNPs. Only pairs with `r2>=0.1` are included. Download [here](https://jh.box.com/s/6hkedzszugebe5ycwh4gt5hjw5sbjdkc). 
* LD score calculated from 1000G Phase 3 European sample, with ~1.2 million HapMap3 SNPs. Only `r2` values greater than 0.01 are included.

Use [PLINK](https://www.cog-genomics.org/plink2) to extract HapMap3 SNPs from 1000G EUR
```
for i in `seq 1 22`;do
        /dcl01/chatterj/data/gqi/plink \
        --bfile ../1000G.EUR.QC.${i} \
        --extract /dcl01/chatterj/data/gqi/1000G_for_simulation/snplist/1000G_EUR_HapMap3_snplist_chr${i}.txt \
        --make-bed \
        --out 1000G.EUR.QC.HM3.${i}
done
```
