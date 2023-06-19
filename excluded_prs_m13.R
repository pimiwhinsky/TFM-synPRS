
system("plink --bfile merged-bfiles --allow-no-sex --remove exclude-A60.txt --make-bed --out mydata-excluded-60")
system("plink --bfile mydata-excluded-60 --allow-no-sex --recode --out mydata-excluded-60")
system("plink --file mydata-excluded-60 --allow-no-sex --pca --out exclude60")
system("plink --bfile mydata-excluded-60 --allow-no-sex --assoc --out assoc-60")

##Age > 60 Ctrls

sp_age <- data.frame(sp_age$IID, sp_age$`Pheno 2 (+ restr)`, sp_age$AGELast, sp_age$AGEFirst)
colnames(sp_age) <- c("IID", "Pheno2", "AgeLast", "AgeFirst")
head(sp_age)

filtro1 <- sp_age$Pheno2 == 1 & sp_age$AgeLast >= 60
df1 <- data.frame(IID = sp_age$IID[filtro1], Pheno = sp_age$Pheno[filtro1], Age = sp_age$AgeLast[filtro1])

filtro2 <- sp_age$Pheno2 == 2 & !is.na(sp_age$AgeFirst)
df2 <- data.frame(IID = sp_age$IID[filtro2], Pheno = sp_age$Pheno[filtro2], Age = sp_age$AgeFirst[filtro2])

filtro3 <- sp_age$Pheno2 == -9 & sp_age$AgeLast
df3 <- data.frame(IID = sp_age$IID[filtro3], Pheno = sp_age$Pheno[filtro3], Age = sp_age$AgeLast[filtro3])
  
df_final <- rbind(df1, df2, df3)
df_final

write.table(df_final, "age60_dataframe.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

IDs_sel <- sp_age$IID[filtro1]
filtro_no_sel <- sp_age$Pheno2 == 1 & sp_age$AgeLast < 60
IDs_no_sel <- sp_age$IID[filtro_no_sel]
write.table(IDs_no_sel, "IDs_no_sel.txt", row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(df_final[, c("IID", "Pheno")], "Pheno2-60.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

##PRS

library(data.table)
covariate <- fread("Pheno2-60.txt")
pcs <- fread("exclude60.eigenvec", header=F)
colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
fwrite(pcs,"synPRS60.covariate", sep="\t")

system("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded-60 --binary-target T --pheno Pheno2-60.txt --cov synPRS60.covariate --stat BETA --no-clump --fastscore --print-snp --out PRSresults-m13")

##COMPARACIÓ PRS X PHENO

prs <- read.table("PRSresults-m13.best", header=T)
AD <- read.table("Pheno2-60.txt", header=T)
head(AD)
dat <- merge(prs, AD, age, by=c("FID","IID"))
dat$Pheno <-as.factor(dat$Pheno)
dat_filtered <- dat[dat$Pheno != -9,]
head(dat_filtered)

###Subset i mean

subset_pheno <- dat[, c("IID", "PRS", "Pheno")]
subset_pheno_filtered <- subset_pheno[subset_pheno$Pheno != -9,]
mean_prs <- aggregate(PRS ~ Pheno, data = subset_pheno_filtered, FUN = mean)
mean_sd <- aggregate(PRS ~ Pheno, data = subset_pheno_filtered, FUN = sd)

print(mean_prs)

age_ctrl_p1 <- pheno1$sp_age.AGELast[pheno1$sp_age.Pheno.1....restr. == 1]
head(age_ctrl_p1)
sd(age_ctrl_p1)
age_AD_p1 <- pheno1$sp_age.AGEFirst[pheno1$sp_age.Pheno.1....restr. == 2]
head(age_AD_p1)
sd(age_AD_p1)

shapiro.test(age_ctrl_p1)
shapiro.test(age_AD_p1)
wilcox.test(age_ctrl_p1, age_AD_p1)
t.test(age_ctrl_p1, age_AD_p1)

age <- read.table("age60_dataframe.txt", header=T)
head(age)
summary_age <- aggregate(Age ~ Pheno, age, summary)
print(summary_age)
mean_sd <- aggregate(Age ~ Pheno, data = age, FUN = sd)
print(mean_sd)

shapiro.test(age$Age[age$Pheno == 1])
shapiro.test(age$Age[age$Pheno == 2])
wilcox.test(age$Age[age$Pheno == 1], age$Age[age$Pheno == 2])
t.test(age_ctrl_p1, age_AD_p1)

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
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  labs(x = "Phenotype", y = "Synaptic PRS Score") + 
  theme_minimal() +
  scale_x_discrete(labels = c("Control", "LOAD"))
ggsave("bpM13.png", width = 5, height = 4, dpi = 300)

##CORRELATION AGE/PRS

age <-read.table("age60_dataframe.txt", header=TRUE)
age <-data.frame(age)
head(age)
df_filtered <- subset(age, Pheno != -9)
head(df_filtered)
df_merged <- merge(df_filtered, df2, by = "IID")
head(df_merged)
df_merged$Age <-as.numeric(df_merged$Age)
df_merged$PRS <-as.numeric(df_merged$PRS)
df_merged$Pheno <-as.factor(df_merged$Pheno)

correlation <- cor(df_merged$PRS, df_merged$Age)
correlation
resultado <- cor.test(df_merged$PRS, df_merged$Age)
print(resultado$estimate)

# Mostrar el valor p
print(resultado$p.value)
print(resultado$estimate)

ggplot(df_merged, aes(x = Age, y = PRS, color = Pheno)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Age", y = "Synaptic PRS") +
  theme_minimal() +
  scale_color_manual(values = c("orange", "cadetblue2"),
                     labels = c("Control", "LOAD")) +
  guides(color = guide_legend(title = " ")) +
  annotate("rect", xmin = 85 - 4, xmax = 85 + 4, ymin = -0.05 - 0.008, ymax = -0.05 + 0.008,
           fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 85, y = -0.05, label = paste("P = 0.019"), size = 3, color = "black")
ggsave("agePRS.png", width = 5, height = 4, dpi = 300)

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
legend("right", legend = paste0("Synaptic PRS-Model 13: AUC = ", auc_value), bty="n")

ggsave("rocM13.png", width = 5, height = 4, dpi = 300)
png("plot.png", width = 800, height = 600, res = 300)
