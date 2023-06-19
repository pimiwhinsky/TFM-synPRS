
system("plink --bfile merged-bfiles --allow-no-sex --remove exclude.txt --make-bed --out mydata-excluded")
system("plink --bfile mydata-excluded --allow-no-sex --recode --out mydata-excluded")
system("plink --bfile mydata-excluded --allow-no-sex --assoc --out assoc-p2")
system("plink --file mydata-excluded --allow-no-sex --pca")

##Age as a cov

library(readxl)
sp_age <- read_excel("sp-age.xlsx")
View(sp_age)

cov.age <- data.frame(IID = sp_age$IID,
                      Pheno = sp_age$`Pheno 2 (+ restr)`,
                      Age = ifelse(sp_age$`Pheno 2 (+ restr)` == 1, sp_age$AGELast, sp_age$AGELast))
View(cov.age)
summary_age <- aggregate(Age ~ Pheno, cov.age, summary)
print(summary_age)
write.table(cov.age, file = "cov_age_sex_p2.txt", row.names = FALSE)

##PRS

library(data.table)
covariate <- fread("cov_age_sex_p2.txt")
head(covariate)
pcs <- fread("plink.eigenvec", header=F)
View(pcs)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
cov <- merge(covariate, pcs)
cov$Age <- as.numeric(cov$Age)
fwrite(cov,"synPRS.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded --binary-target T --pheno pheno2.txt --cov synPRS.covariate --stat BETA --fastscore --print-snp --out PRSresults_m5")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults_m5.best", header=T)
covar <- read.table("synPRS.covariate", header=T)
dat <- merge(prs, covar, by=c("FID","IID"))
View(dat)
dat$Pheno <-as.factor(dat$Pheno)
dat$Sex <-as.factor(dat$Sex)
dat$Age <-as.numeric(dat$Age)
dat_filtered <- dat[c(dat$Pheno != -9, dat$Sex != 0),]
head(dat_filtered)

###Subset i mean

subset_pheno <- dat_filtered[, c("IID", "PRS", "Pheno", "Sex", "Age")]
mean_prs <- aggregate(PRS ~ Pheno, data = subset_pheno, FUN = mean)
print(mean_prs)
mean_prs_sex <- aggregate(PRS ~ Pheno + Sex, data = subset_pheno_filtered, FUN = mean)
print(mean_prs_sex)

summary_age <- aggregate(Age ~ Pheno, subset_pheno_filtered, summary)
print(summary_age)
summary_age_sex <- aggregate(Age ~ Sex, subset_pheno_filtered, summary)
print(summary_age_sex)
media_edad <- aggregate(Age ~ Pheno + Sex, data = subset_pheno_filtered, FUN = mean)
print(media_edad)
mean(subset_pheno_filtered$Age)

summary_sex <- tapply(subset_pheno_filtered$Pheno, subset_pheno_filtered$Sex, summary)
print(summary_sex)
summary_sex_nf <- tapply(dat_filtered$Pheno, dat_filtered$Sex, summary)
print(summary_sex_nf)

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
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype - Model 5") +
  scale_x_discrete(labels = c("Control", "LOAD"))

ggplot(dat_filtered2, aes(x = Pheno, y = PRS, fill = factor(Sex))) + 
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(alpha = 0.5, position = position_jitterdodge(), size = 3, shape = 4, alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"), name = "Sex", labels = c("M", "F")) +
  labs(x = "Phenotype", y = "Synaptic PRS Score", fill = "Sex")  +
  ggtitle("Boxplot of Synpatic PRS Score by Phenotype and Sex - Model 5") +
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
legend("right", legend = paste0("Synaptic PRS-Model 5: AUC = ", auc_value), bty="n")