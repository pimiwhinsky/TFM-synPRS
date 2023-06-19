
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


##MISSING DATA, MAFS, ASSOC, HWE, EIGENV

options <- c("F:\\dades\\nas\\GRACE\\TopMED_forEADB\\DEGESCO_inEADB\\SPIN\\coverted-bfiles\\PRS\\excluded_PRS\\synPRS\\pheno1\\", "F:\\dades\\nas\\GRACE\\TopMED_forEADB\\DEGESCO_inEADB\\SPIN\\coverted-bfiles\\PRS\\excluded_PRS\\synPRS\\PHENO2\\")
Eigenvalue_docnames <- c("plink1.eigenval","plink2.eigenval")
Eigenvector_docnames <- c("plink1.eigenvec","plink2.eigenvec")
c = 1

for (i in options){
  system(paste0("plink --bfile ",i,"mydata-excluded --allow-no-sex --missing --out ",i,"callrate_bfiles_p",str(c)))
  system(paste0("plink --bfile ",i,"mydata-excluded --allow-no-sex --freq case-control --out ",i,"mafs_p",str(c)))
  system(paste0("plink --bfile ",i,"mydata-excluded --allow-no-sex --assoc --out ",i,"assoc_p",str(c)))
  system(paste0("plink --bfile ",i,"mydata-excluded --allow-no-sex --hardy --out ",i,"hwe_bfiles_p",str(c)))
  system(paste0("plink --file mydata-excluded --allow-no-sex --pca"))
  system("mv plink.eigenvec, ",Eigenvalue_docnames[c])
  system("mv plink.eigenval, ", Eigenvector_docnames[c])
  
  c = c+1
}


# PRS2.0

options <- c("plink1.eigenvec","plink2.eigenvec")
docs <- c("synPRS.covariate","synPRS-nocov.covariate")
phenos <- c("pheno1.txt","pheno2.txt")

c = 0
d = 0
for (file in options){
  d = d+1
  pcs <- fread(file, header=F)
  colnames(pcs) <- c("FID","IID", paste0("PC",1:8))
  for (i in range(length(docs))){
    if (i == 1){
      covariate <- fread("cov.txt")
      cov <- merge(covariate, pcs)
      fwrite(cov,docs[1], sep="\t")
      for (j in 1:2){
        c=c+1
        if (j == 1){
          system(paste0("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded",str(d)," --binary-target T --pheno ",phenos[j]," --cov ",docs[i]," --stat BETA --fastscore --print-snp --out PRSresults-m",str(c)))
          
        }else{
          system(paste0("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded",str(d)," --binary-target T --pheno ",phenos[j]," --cov ",docs[i]," --stat BETA --no-clump --fastscore --print-snp --out PRSresults-m",str(c)))
        }
      }
      
    }else{
      fwrite(pcs,docs[2], sep="\t")
      for (j in 1:2){
        c=c+1
        if (j == 1){
          system(paste0("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded",str(d)," --binary-target T --pheno ",phenos[j]," --cov ",docs[i]," --stat BETA --fastscore --print-snp --out PRSresults-m",str(c)))
        }else{
          system(paste0("Rscript PRSice.R --prsice PRSice_win64.exe --base summstats.txt --target mydata-excluded",str(d)," --binary-target T --pheno ",phenos[j]," --cov ",docs[i]," --stat BETA --no-clump --fastscore --print-snp --out PRSresults-m",str(c)))
        }
      }
      
    }
    
  }
}

###IDs

install.packages("writexl")
excel_data <- read_excel("PhenoEADB-160223.xlsx", sheet = 1)
head(excel_data)
filtered_data <- subset(excel_data, Pheno1...2 == 1, select = c("IDNUM", "SP...3"))
write.csv(filtered_data, file="pheno2ID.csv", row.names=FALSE)

##ASSOCIATIONS

library(readr)
library(ggplot2)
library(dplyr)

pheno2ID_datos <- read_delim("pheno2ID-datos.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(pheno2ID_datos)

df2 <-read.table("PRSresults-m13.best", header = TRUE)
head(df2)
df3 <-read.table("pheno2.txt", header = TRUE)
head(df3)
df4 <-merge(df2, df3, by=c("FID", "IID"))
head(df4)
df5 <- subset(df4, Pheno == 1)

df <-data.frame(pheno2ID_datos$IDNUM, pheno2ID_datos$Months, pheno2ID_datos$MMSEFirst, pheno2ID_datos$FCSRTTOTALNORM, pheno2ID_datos$Ratio4240, pheno2ID_datos$CSFPTAU_CSF_E0300, pheno2ID_datos$CSFTAU_CSF_E0200, df5$PRS)
colnames(df) <-c("ID", "Months", "MMSE", "FCSRT", "Ratio", "CSF300", "CSF200", "PRS")
head(df)

#FCSRT

df6 <-na.omit(data.frame(df$ID, df$Months, df$FCSRT, df$PRS))
View(df6)
summary(df6)

model <- lm(df.FCSRT ~ df.PRS + df.Months, data = df6)
summary(model)

tertiles <- quantile(df6$df.PRS, probs = c(0, 1/3, 2/3, 1))
df6 <- df6 %>% mutate(Tertiles = cut(df.PRS, breaks = tertiles, labels = c("T1", "T2", "T3")))
head(df6)
df6_subset <- na.omit(df6 %>% select(Tertiles, df.FCSRT))

rango_tertiles <- paste(round(tertiles[-length(tertiles)], 2), "-", round(tertiles[-1], 2))

ggplot(df6_subset, aes(x = Tertiles, y = df.FCSRT)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  
  scale_x_discrete(labels = rango_tertiles) +  
  xlab("LOAD PRS Model 13") +
  ylab("FCSRT") +
  theme_minimal() +
  annotate("rect", xmin = 1 - 0.4, xmax = 1 + 0.4, ymin = 16 - 0.5, ymax = 16 + 0.4,
           fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 1, y = 16, label = paste("p-v = 0.2628"), size = 4, color = "black")
ggsave("FCSRT.png", width = 5, height = 4, dpi = 300)

#MMSE

df7 <-na.omit(data.frame(df$ID, df$Months, df$MMSE, df$PRS))
head(df7)
summary(df7)

model <- lm(df.MMSE ~ df.PRS + df.Months, data = df7)
summary(model)

tertiles <- quantile(df7$df.PRS, probs = c(0, 1/3, 2/3, 1))
df7 <- df7 %>% mutate(Tertiles = cut(df.PRS, breaks = tertiles, labels = c("T1", "T2", "T3")))
View(df7)
df7_subset <- na.omit(df7 %>% select(Tertiles, df.MMSE))

rango_tertiles <- paste(round(tertiles[-length(tertiles)], 2), "-", round(tertiles[-1], 2))

ggplot(df7_subset, aes(x = Tertiles, y = df.MMSE)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  
  scale_x_discrete(labels = rango_tertiles) +  
  xlab("LOAD PRS Model 13") +
  ylab("MMSE") +
  theme_minimal() +
  annotate("rect", xmin = 3.1 - 0.4, xmax = 3.1 + 0.4, ymin = 5 - 0.7, ymax = 5 + 0.6,
           fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 3.1, y = 5, label = paste("p-v = 0.8777"), size = 4, color = "black")
ggsave("MMSE.png", width = 5, height = 4, dpi = 300)

#AB

df8 <-na.omit(data.frame(df$ID, df$Months, df$Ratio, df$PRS))
df8$df.Ratio <- gsub(",", ".", df8$df.Ratio)
df8$df.Ratio <- as.numeric(df8$df.Ratio)
head(df8)

model <- lm(df.Ratio ~ df.PRS + df.Months, data = df8)
summary(model)

tertiles <- quantile(df8$df.PRS, probs = c(0, 1/3, 2/3, 1))
df8 <- df8 %>% mutate(Tertiles = cut(df.PRS, breaks = tertiles, labels = c("T1", "T2", "T3")))
View(df8)
df8_subset <- na.omit(df8 %>% select(Tertiles, df.Ratio))

rango_tertiles <- paste(round(tertiles[-length(tertiles)], 2), "-", round(tertiles[-1], 2))

ggplot(df8_subset, aes(x = Tertiles, y = df.Ratio)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  
  scale_x_discrete(labels = rango_tertiles) +  
  xlab("LOAD PRS Model 13") +
  ylab("Aβ42/40") +
  theme_minimal() +
  annotate("rect", xmin = 0.5 - 0.4, xmax = 0.5 + 0.4, ymin = 0.065 - 0.001, ymax = 0.065 + 0.0009,
         fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 0.5, y = 0.065, label = paste("p-v = 0.6091"), size = 4, color = "black")
ggsave("AB.png", width = 5, height = 4, dpi = 300)

#pTAU

df9 <-na.omit(data.frame(df$ID, df$Months, df$CSF300, df$PRS))
head(df9)

model <- lm(df.CSF300 ~ df.PRS + df.Months, data = df9)
summary(model)

tertiles <- quantile(df9$df.PRS, probs = c(0, 1/3, 2/3, 1))
df9 <- df9 %>% mutate(Tertiles = cut(df.PRS, breaks = tertiles, labels = c("T1", "T2", "T3")))
df9_subset <- na.omit(df9 %>% select(Tertiles, df.CSF300))

rango_tertiles <- paste(round(tertiles[-length(tertiles)], 2), "-", round(tertiles[-1], 2))

ggplot(df9_subset, aes(x = Tertiles, y = df.CSF300)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  
  scale_x_discrete(labels = rango_tertiles) +  
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1100, by = 100)) + 
  xlab("LOAD PRS Model 13") +
  ylab("pTau") +
  theme_minimal() +
  annotate("rect", xmin = 1 - 0.4, xmax = 1 + 0.4, ymin = 1100 - 10, ymax = 1100 + 10,
           fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 1, y = 1100, label = paste("p-v = 0.2129"), size = 4, color = "black")
ggsave("tau300.png", width = 5, height = 4, dpi = 300)

#TAU

df10 <-na.omit(data.frame(df$ID, df$Months, df$CSF200, df$PRS))
head(df10)

model <- lm(df.CSF200 ~ df.PRS + df.Months, data = df10)
summary(model)

tertiles <- quantile(df10$df.PRS, probs = c(0, 1/3, 2/3, 1))
df10 <- df10 %>% mutate(Tertiles = cut(df.PRS, breaks = tertiles, labels = c("T1", "T2", "T3")))
df10_subset <- na.omit(df10 %>% select(Tertiles, df.CSF200))

rango_tertiles <- paste(round(tertiles[-length(tertiles)], 2), "-", round(tertiles[-1], 2))

ggplot(df10_subset, aes(x = Tertiles, y = df.CSF200)) +
  geom_jitter(position = position_jitter(width = 0.2), size = 3, shape = 4, alpha = 1) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +  
  scale_x_discrete(labels = rango_tertiles) +  
  scale_y_continuous(limits = c(0, 1100), breaks = seq(0, 1100, by = 100)) + 
  xlab("LOAD PRS Model 13") +
  ylab("tTau") +
  theme_minimal() +
  annotate("rect", xmin = 3 - 0.4, xmax = 3 + 0.4, ymin = 50 - 35, ymax = 50 + 25,
            fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 3, y = 50, label = paste("p-v = 0.8981"), size = 4, color = "black")
ggsave("tau200.png", width = 5, height = 4, dpi = 300)

#Age

age <-read.table("cov_age_p2.txt", header=TRUE)
age <-data.frame(age)
head(age)
df_filtered <- subset(age, Pheno != -9)
head(df_filtered)
df_merged <- merge(df_filtered, df2, by = "IID")
head(df_merged)
df_merged$Age <-as.numeric(df_merged$Age)
df_merged$PRS <-as.numeric(df_merged$PRS)
df_merged$Pheno <-as.factor(df_merged$Pheno)

aggregate(Age ~ Pheno, data = df_merged, FUN = mean)
aggregate(PRS ~ Pheno, data = df_merged, FUN = mean)

library(ggplot2)
library(dplyr)
correlation <- cor(df_merged$PRS, df_merged$Age)
correlation

ggplot(df_merged, aes(x = Age, y = PRS, color = Pheno)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  xlab( "Age") +
  ylab("Synaptic PRS Model 8") +
  theme_minimal() +
  scale_color_manual(values = c("orange", "cadetblue2"),
                     labels = c("Control", "LOAD")) +
  guides(color = guide_legend(title = " ")) +
  annotate("rect", xmin = 55 - 7, xmax = 55 + 7, ymin = 0.2 - 0.008, ymax = 0.2 + 0.008,
         fill = "white", color = "red", linetype = "solid") +
  annotate("text", x = 55, y = 0.2, label = paste("Corr = -0.0267"), size = 4, color = "black")
ggsave("age.png", width = 5, height = 4, dpi = 300)

