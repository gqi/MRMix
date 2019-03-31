# This script implements LD clumping with threshold r2<r2thr.
# @param sumstats Data frame containing chr, rsid, bp, z, pval
# @param LDinfo r2 coefficient (square of LD coefficient) of a SNP set including those in sumstats
# @param r2thr r2 threshold for LD clumping
ld_clump_MRsim = function(sumstats, LDinfo, pthr = 5e-8, r2thr = 0.1){

    tempdata = sumstats %>% filter(pval<pthr) %>% rename(CHR = chr, SNP = rsid, BP = bp, P = pval)

    snps_clumped_all = data.frame(CHR = vector(length=0), SNP = vector(length=0),
                                  BP = vector(length=0), z = vector(length=0), P = vector(length=0))
    if (r2thr>0.1){
        for (chrnum in 1:22){
            LDinfo[[chrnum]] = LDinfo[[chrnum]] %>% filter(R2>r2thr)
        }
    }

    for (chrnum in 1:22){
        print(chrnum)
        if (sum(tempdata$CHR==chrnum)>0){
            tempdata_chr = tempdata %>% filter(CHR==chrnum)
            if (nrow(tempdata_chr)>1){
                tempdata_chr = tempdata_chr[order(tempdata_chr$P,decreasing=FALSE),]
                temp_pruned = NULL
                while (nrow(tempdata_chr)>0){
                    temp_pruned = rbind(temp_pruned,tempdata_chr[1,])
                    if (nrow(tempdata_chr)>1){
                        ind = (tempdata_chr$SNP %in% c(LDinfo[[chrnum]]$SNP_A[LDinfo[[chrnum]]$SNP_B==tempdata_chr$SNP[1]],LDinfo[[chrnum]]$SNP_B[LDinfo[[chrnum]]$SNP_A==tempdata_chr$SNP[1]]))
                        tempdata_chr = tempdata_chr[-c(1,which(ind)),]
                    } else{
                        tempdata_chr = tempdata_chr[-1,]
                    }
                }
                snps_clumped_all = rbind(snps_clumped_all, temp_pruned)
            } else{
                snps_clumped_all = rbind(snps_clumped_all, tempdata_chr)
            }
        }
    }

    return(snps_clumped_all)
}
