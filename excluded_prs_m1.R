
system("plink --bfile bfile1 --bmerge bfile2.bed bfile2.bim bfile2.fam --allow-no-sex --make-bed --out mergebfile1")
system("plink --bfile mergebfile1 --bmerge bfile3.bed bfile3.bim bfile3.fam --allow-no-sex --make-bed --out mergebfile2")
system("plink --bfile mergebfile2 --bmerge bfile4.bed bfile4.bim bfile4.fam --allow-no-sex --make-bed --out mergebfile3")
system("plink --bfile mergebfile3 --bmerge bfile5.bed bfile5.bim bfile5.fam --allow-no-sex --make-bed --out merged-bfiles")

system("plink --bfile merged-bfiles --allow-no-sex --remove exclude-A.txt --make-bed --out mydata-excluded-A")
system("plink --bfile mydata-excluded-A --allow-no-sex --recode --out mydata-excluded-A")
system("plink --bfile mydata-excluded-A --allow-no-sex --assoc --out assoc-A-p1")

##CALL RATIO

system("plink --bfile mydata-excluded-A --allow-no-sex --missing --out callrate-bfiles")

##MAF

system("plink --bfile mydata-excluded-A --allow-no-sex --freq case-control --out mafs_p1")

##HWE

system("plink --bfile mydata-excluded-A --allow-no-sex --hardy --out hwe_bfiles_p1")

##EIGENVECTORS

system("plink --file mydata-excluded-A --allow-no-sex --pca")

##PRS

library(data.table)
covariate <- fread("cov_excluded-A.txt")
pcs <- fread("plink.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
cov <- merge(covariate, pcs)
fwrite(cov,"synPRS-A.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded-A --binary-target T --pheno pheno1-A.txt --cov synPRS-A.covariate --stat BETA --fastscore --print-snp --out PRSresults-m1A")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m1A.best", header=T)
AD <- read.table("pheno1-A.txt", header=T)
covar <- read.table("synPRS-A.covariate", header=T)
covar$Sex <- as.factor(covar$Sex)
dat <- merge(merge(prs, AD), covar, by=c("FID","IID"))
dat$Pheno <-as.factor(dat$Pheno)
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

summary_age <- aggregate(Age ~ Pheno, subset_pheno_filtered, summary)
print(summary_age)

summary_age_sex <- aggregate(Age ~ Sex, subset_pheno_filtered, summary)
print(summary_age_sex)

media_edad <- aggregate(Age ~ Pheno + Sex, data = subset_pheno_filtered, FUN = quantile)
print(media_edad)

summary_sex <- tapply(subset_pheno_filtered$Pheno, subset_pheno_filtered$Sex, summary)
print(summary_sex)

summary_sex_nf <- tapply(dat$Pheno, dat$Sex, summary)
print(summary_sex_nf)

mean(subset_pheno_filtered$Age)

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
  labs(x = "Phenotype", y = "Synaptic PRS Score") + 
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype - Model 1") +
  scale_x_discrete(labels = c("Control", "LOAD"))

ggplot(dat_filtered2, aes(x = Pheno, y = PRS, fill = factor(Sex))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(), size = 3, shape = 4, alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Sex", labels = c("M", "F")) +
  labs(x = "Phenotype", y = "Synaptic PRS Score", fill = "Sex")  +
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype and Sex - Model 1") +
  scale_x_discrete(labels = c("Control", "LOAD"))

##CORRELATION AGE/PRS

correlation <- cor(subset_pheno_filtered$PRS, subset_pheno_filtered$Age)
correlation

ggplot(subset_pheno_filtered, aes(x = Age, y = PRS, color = Pheno)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  geom_text(aes(label = paste("Corr:", round(correlation, 4))),
            x = Inf, y = -Inf, hjust = 1, vjust = 0,
            size = 4, color = "black") +
  labs(x = "Age", y = "Synaptic PRS") +
  theme_minimal() +
  scale_color_manual(values = c("orange", "cadetblue2"),
                     labels = c("Control", "LOAD")) +
  guides(color = guide_legend(title = " "))

###pROC

subset_pheno_filtered$Pheno <- relevel(factor(subset_pheno_filtered$Pheno), ref = "1")

library(pROC)
roc_obj <- roc(subset_pheno_filtered$Pheno, subset_pheno_filtered$PRS, percent=TRUE)
auc_value <- round(auc(roc_obj), 2)
auc_value
ci <- ci.auc(roc_obj)
ci
plot(roc_obj, print.auc=FALSE, auc.polygon=TRUE, print.thres=FALSE, grid=c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9))
legend("right", legend = paste0("Synaptic PRS-Model 1: AUC = ", auc_value), bty="n")