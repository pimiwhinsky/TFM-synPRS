
system("plink --vcf SPIN_chr1_imputed_GSA.vcf.gz --make-bed --out chr19")
system("plink --bfile chr19 --allow-no-sex --extract rsAPOE.txt --make-bed --out chr19-rsAPOE")

system("plink --vcf SPIN_chr2_imputed_GSA.vcf.gz --make-bed --out chr2")
system("plink --bfile chr2 --allow-no-sex --extract rsBIN1.txt --make-bed --out chr2-rsBIN1")

system("plink --vcf SPIN_chr8_imputed_GSA.vcf.gz --make-bed --out chr8")
system("plink --bfile chr8 --allow-no-sex --extract rsPTK2B.txt --make-bed --out chr8-rsPTK2B")

system("plink --vcf SPIN_chr11_imputed_GSA.vcf.gz --make-bed --out chr11")
system("plink --bfile chr11 --allow-no-sex --extract rschr11.txt --make-bed --out chr11-rschr11")

system("plink --vcf SPIN_chr17_imputed_GSA.vcf.gz --make-bed --out chr17")
system("plink --bfile chr17 --allow-no-sex --extract rsMINK1.txt --make-bed --out chr17-rsMINK1")

##MERGE

system("plink --bfile bfile1 --bmerge bfile2.bed bfile2.bim bfile2.fam --allow-no-sex --make-bed --out mergebfile1")
system("plink --bfile mergebfile1 --bmerge bfile3.bed bfile3.bim bfile3.fam --allow-no-sex --make-bed --out mergebfile2")
system("plink --bfile mergebfile2 --bmerge bfile4.bed bfile4.bim bfile4.fam --allow-no-sex --make-bed --out mergebfile3")
system("plink --bfile mergebfile3 --bmerge bfile5.bed bfile5.bim bfile5.fam --allow-no-sex --make-bed --out merged-bfiles")

system("plink --bfile merged-bfiles --allow-no-sex --remove exclude.txt --make-bed --out mydata-excluded")

system("plink --bfile mydata-excluded --allow-no-sex --recode --out mydata-excluded")

##CALL RATIO

system("plink --bfile mydata-excluded --allow-no-sex --missing --out callrate_bfiles_p2")
  
##MAF

system("plink --bfile mydata-excluded --allow-no-sex --freq case-control --out mafs_p2")

##ASSOC

system("plink --bfile mydata-excluded --allow-no-sex --assoc --out assoc_p2")

##HWE

system("plink --bfile mydata-excluded --allow-no-sex --hardy --out hwe_bfiles_p2")

##EIGENVECTORS

system("plink --file mydata-excluded --allow-no-sex --pca")

##PRS

library(data.table)
pcs <- fread("plink.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
fwrite(pcs,"synPRS-nocov.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded --binary-target T --pheno pheno2.txt --cov synPRS-nocov.covariate --stat BETA --perm 10000 --no-clump --fastscore --print-snp --out PRSresults-m8")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m8.best", header=T)
AD <- read.table("pheno2.txt", header=T)
dat <- merge(prs, AD, by=c("FID","IID"))
dat$Pheno <-as.factor(dat$Pheno)
head(dat)

##Sex analysis

sexp2 <- read.table("cov_age_sex_p2.txt", header=T)
head(sexp2)
View(sexp2)
sexp2$Sex <-as.factor(sexp2$Sex)
table(sexp2$Sex)
datos_subset2 <- subset(sexp2, Pheno %in% c(1, 2))

tabla_contingencia <- table(datos_subset2$Pheno, datos_subset2$Sex)
print(tabla_contingencia)

test_chi <- chisq.test(tabla_contingencia)
print(test_chi)

sex <- read.table("cov_excluded-A.txt", header=T)
AD <-read.table("pheno1-A.txt", header=T)
sexp1 <- merge(sex, AD, by=c("FID","IID"))
head(sexp1)
sexp2$Sex <-as.factor(sexp1$Sex)
datos_subset <- subset(sexp1, Pheno %in% c(1, 2))

tabla_contingencia <- table(datos_subset$Pheno, datos_subset$Sex)
print(tabla_contingencia)

test_chi <- chisq.test(tabla_contingencia)
print(test_chi)

###Subset i mean

subset_pheno <- dat[, c("IID", "PRS", "Pheno")]
subset_pheno_filtered <- subset_pheno[subset_pheno$Pheno != -9,]
mean_prs <- aggregate(PRS ~ Pheno, data = subset_pheno_filtered, FUN = mean)
print(mean_prs)

###Quantiles

prs_cases <- subset_pheno_filtered$PRS[subset_pheno_filtered$Pheno == 2]
prs_controls <- subset_pheno_filtered$PRS[subset_pheno_filtered$Pheno == 1]
quantiles_cases <- quantile(prs_cases, probs = seq(0, 1, by = 0.25))
quantiles_controls <- quantile(prs_controls, probs = seq(0, 1, by = 0.25))
print(quantiles_cases)
print(quantiles_controls)

###Levene i ttest/wilcox

library(car)
shapiro_test1 <- shapiro.test(prs_cases)
print(shapiro_test1)
shapiro_test2 <- shapiro.test(prs_controls)
print(shapiro_test2)
levene_test <- leveneTest(PRS ~ Pheno, data = subset_pheno_filtered)
print(levene_test)
t_test <- t.test(PRS ~ Pheno, data = subset_pheno_filtered)
print(t_test)
wilcox.test(PRS ~ Pheno, data = subset_pheno_filtered)

###Boxplot Pheno/PRS

library(ggplot2)
ggplot(subset_pheno_filtered, aes(x = factor(Pheno), y = PRS)) + 
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 0.5) +
  labs(x = "Phenotype", y = "Synaptic PRS Score") + 
  theme_minimal() +
  scale_x_discrete(labels = c("Control", "LOAD"))
ggsave("bpM8.png", width = 5, height = 4, dpi = 300)

###pROC

subset_pheno_filtered$Pheno <- relevel(factor(subset_pheno_filtered$Pheno), ref = "1")
levels(subset_pheno_filtered$Pheno)

library(pROC)
roc_obj <- roc(subset_pheno_filtered$Pheno, subset_pheno_filtered$PRS, percent=TRUE)
auc_value <- round(auc(roc_obj), 2)
auc_value
ci <- ci.auc(roc_obj)
ci
plot(roc_obj, print.auc=FALSE, auc.polygon=TRUE, print.thres=FALSE, grid=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
legend("right", legend = paste0("Synaptic PRS-Model 8: AUC = ", auc_value), bty="n")

