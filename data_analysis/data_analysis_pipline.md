# MRMix data analysis pipeline

To reproduce the results in our paper (Qi and Chatterjee, Nat Comms 2018), please download the following files:

* A PLINK `bim` file for HapMap 3 SNPs with MAF>0.05 1000 Genomes Phase 3 European sample (around 1.07 million SNPs). Download [here](https://www.dropbox.com/s/rga3na4zdq9oaxq/1000G.EUR.QC.HM3.MERGE.maf05.bim?dl=0).
* Linkage disequilibrium (LD) coefficients calculated from 1000G Phase 3 European sample, with ~1.2 million HapMap3 SNPs. Only pairs with `r2>=0.1` are included. Download [here](https://www.dropbox.com/s/yby3igxx6wwy120/LDinfo_with_R2.rda?dl=0). 
* LD score calculated from 1000G Phase 3 European sample, with ~1.2 million HapMap3 SNPs. Only `r2` values greater than 0.01 are included. Download [here](https://www.dropbox.com/s/454l3q1ll9wpwvd/ldscore_thr_r0.1.rda?dl=0).

Code for preprocessing and links to summary level data are provided [here](https://github.com/gqi/MRMix/blob/master/data_analysis/preprocess_summary_stats.R). 

Code for Mendelian randomziation analysis using `MRMix` and `MendelianRandomization` package is provided [here](https://github.com/gqi/MRMix/blob/master/data_analysis/MR_realdata_analysis.R).

In our paper we used [our own code](https://github.com/gqi/MRMix/blob/master/data_analysis/ld_clump_MRsim.R) for LD clumping with `r2<=0.1`. [PLINK](https://www.cog-genomics.org/plink2) can also be used for LD clmping. 

Note that the body mass index (BMI) summary statistics were updated on 6/25/2018. Our results that involves BMI are based on an earlier version of the data and may not be exactly reproduced. 
