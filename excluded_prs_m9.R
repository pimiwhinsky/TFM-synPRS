
system("plink --vcf SPIN_chr1_imputed_GSA.vcf.gz --make-bed --out chr19")
system("plink --bfile chr19 --allow-no-sex --extract rsAPOE.txt --make-bed --out bfile5")

system("plink --bfile bfile5 --allow-no-sex --remove exclude.txt --make-bed --out mydata-excluded-apoe")

system("plink --bfile mydata-excluded-apoe --allow-no-sex --recode --out mydata-excluded-apoe")

##CALL RATIO

system("plink --bfile mydata-excluded-apoe --allow-no-sex --missing --out callrate_apoe")
  
##MAF

system("plink --bfile mydata-excluded-apoe --allow-no-sex --freq case-control --out mafs_apoe_p1")

##ASSOC

system("plink --bfile mydata-excluded-apoe --allow-no-sex --assoc --out assoc_apoe_p1")

##HWE

system("plink --bfile mydata-excluded-apoe --allow-no-sex --hardy --out hwe_bfile_apoe_p1")

##EIGENVECTORS

system("plink --file mydata-excluded-apoe --allow-no-sex --pca")

##PRS

library(data.table)
covariate <- fread("cov_excluded.txt")
pcs <- fread("plink.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:2))
cov <- merge(covariate, pcs)
fwrite(cov,"syn-apoe.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstat-apoe.txt --target mydata-excluded-A --binary-target T --pheno pheno1-A.txt --cov synPRS-A.covariate --stat BETA --fastscore --print-snp --out PRSresults-m9")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m9.best", header=T)
AD <- read.table("pheno1-A.txt", header=T)
covar <- read.table("synPRS-A.covariate", header=T)
dat <- merge(merge(prs, AD), covar, by=c("FID","IID"))
dat$Pheno <-as.factor(dat$Pheno)
dat$Sex <- as.factor(dat$Sex)
dat_filtered <- dat[dat$Pheno != -9,]
dat_filtered2 <- dat_filtered[dat_filtered$Sex != 0,]
head(dat_filtered2)

###Subset i mean


subset_pheno <- dat_filtered2[, c("IID", "PRS", "Pheno", "Sex", "Age")]
subset_pheno_filtered <- subset_pheno[subset_pheno$Pheno != -9,]
mean_prs <- aggregate(PRS ~ Pheno, data = subset_pheno_filtered, FUN = mean)
print(mean_prs)
mean_prs_sex <- aggregate(PRS ~ Pheno + Sex, data = subset_pheno_filtered, FUN = mean)
print(mean_prs_sex)

###Quantiles

prs_cases <- subset_pheno_filtered$PRS[subset_pheno_filtered$Pheno == 2]
prs_controls <- subset_pheno_filtered$PRS[subset_pheno_filtered$Pheno == 1]
quantiles_cases <- quantile(prs_cases, probs = seq(0, 1, by = 0.25))
quantiles_controls <- quantile(prs_controls, probs = seq(0, 1, by = 0.25))
print(quantiles_cases)
print(quantiles_controls)
quantiles_by_pheno_sex <- aggregate(PRS ~ Pheno + Sex, data = subset_pheno_filtered, 
                                    FUN = function(x) quantile(x, probs = seq(0, 1, by = 0.25)))
print(quantiles_by_pheno_sex)

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
  labs(x = "Phenotype", y = "APOE PRS Score") + 
  ggtitle("Boxplot of APOE PRS Score by Phenotype - Model 9") +
  scale_x_discrete(labels = c("Control", "LOAD"))

ggplot(subset_pheno_filtered, aes(x = Pheno, y = PRS, fill = factor(Sex))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(), size = 3, shape = 4, alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Sex", labels = c("M", "F")) +
  labs(x = "Phenotype", y = "APOE PRS Score", fill = "Sex")  +
  ggtitle("Boxplot of APOE PRS Score by Phenotype and Sex - Model 9") +
  scale_x_discrete(labels = c("Control", "LOAD"))

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
legend("right", legend = paste0("APOE PRS-Model 9: AUC = ", auc_value), bty="n")
