##PRS

library(data.table)
covariate <- fread("cov_excluded.txt")
pcs <- fread("plink.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
cov <- merge(covariate, pcs)
fwrite(cov,"synPRS-A.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded-A --binary-target T --pheno pheno1-A.txt --cov synPRS-A.covariate --stat BETA --no-clump --fastscore --print-snp --out PRSresults-m2A")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m2A.best", header=T)
AD <- read.table("pheno1-A.txt", header=T)
covar <- read.table("synPRS-A.covariate", header=T)
covar$Sex <- as.factor(covar$Sex)
dat <- merge(merge(prs, AD), covar, by=c("FID","IID"))
summary(dat)
dat$Pheno <-as.factor(dat$Pheno)
dat_filtered <- dat[dat$Pheno != -9,]
dat_filtered2 <- dat_filtered[dat_filtered$Sex != 0,]
head(dat_filtered2)

###Subset i mean

subset_pheno <- dat_filtered2[, c("IID", "PRS", "Pheno", "Sex", "Age")]
subset_pheno_filtered <- subset_pheno[subset_pheno$Pheno != -9,]
head(subset_pheno_filtered)
mean_prs <- aggregate(PRS ~ Pheno, data = subset_pheno_filtered, FUN = mean)
print(mean_prs)
mean_prs_sex <- aggregate(PRS ~ Pheno + Sex, data = subset_pheno_filtered, FUN = mean)
print(mean_prs_sex)

media_edad <- aggregate(Age ~ Pheno, data = subset_pheno_filtered, FUN = quantile)
print(media_edad)
media_edad_sexo <- aggregate(Age ~ Pheno + Sex, data = subset_pheno_filtered, FUN = quantile)
print(media_edad_sexo)

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
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype - Model 2") +
  scale_x_discrete(labels = c("Control", "LOAD"))

ggplot(dat_filtered2, aes(x = Pheno, y = PRS, fill = factor(Sex))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(), size = 3, shape = 4, alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Sex", labels = c("M", "F")) +
  labs(x = "Phenotype", y = "Synaptic PRS Score", fill = "Sex")  +
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype and Sex - Model 2") +
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
legend("right", legend = paste0("Synaptic PRS-Model 2: AUC = ", auc_value), bty="n")