##PRS

library(data.table)
pcs <- fread("plink.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
fwrite(pcs,"synPRS-nocov.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstat-apoe.txt --target mydata-excluded-A --binary-target T --pheno pheno1-A.txt --cov synPRS-nocov.covariate --stat BETA --fastscore --print-snp --out PRSresults-m10")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m10.best", header=T)
AD <- read.table("pheno1-A.txt", header=T)
dat <- merge(prs, AD, by=c("FID","IID"))
dat$Pheno <-as.factor(dat$Pheno)
head(dat)

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
  labs(x = "Phenotype", y = "APOE PRS Score") + 
  ggtitle("Boxplot of APOE PRS Score by Phenotype - Model 10") +
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
legend("right", legend = paste0("APOE PRS-Model 10: AUC = ", auc_value), bty="n")
