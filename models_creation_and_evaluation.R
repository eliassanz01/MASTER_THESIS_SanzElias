library(rstatix)
library(car)
library(broom)
library(dplyr)
library(tidyr)
library(stringr)
library(readODS)

# ------------------------------
# 0. DF
# ------------------------------
setwd("/home/admin/Desktop/ELIAS/pafip_PRSCS/")
pafip <- read.csv("pheno.csv", sep="\t")

pafip$p2diag6m <- as.factor(pafip$p2diag6m)
pafip$p3tto3m <- as.factor(pafip$p3tto3m)
pafip$p3tto0 <- as.factor(pafip$p3tto0)
pafip$Alcohol <- as.factor(pafip$Alcohol)
pafip$Cannabis <- as.factor(pafip$Cannabis)
pafip$Tobacco <- as.factor(pafip$Tobacco)

# ------------------------------
# 1. Descriptives
# ------------------------------
summary(pafip)
colnames(pafip)

# Long format
pafip_long <- pafip %>%
  select(BMI0, BMI3M, BMI1a, BMI3a) %>%
  pivot_longer(cols = everything(),
               names_to = "Timepoint",
               values_to = "BMI")

pafip_long$Timepoint <- factor(pafip_long$Timepoint,
                               levels = c("BMI0", "BMI3M", "BMI1a", "BMI3a"),
                               labels = c("Baseline", "3 months", "1 year", "3 years"))

# Create the boxplots
ggplot(pafip_long, aes(x = Timepoint, y = BMI)) +
  geom_boxplot(fill = "#69b3a2", alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.3, color = "gray40") +  # Opcional: puntos individuales
  labs(title = "BMI Evolution in patients experiencing a FEP",
       x = "Timepoint",
       y = "Body Mass Index") +
  theme_minimal(base_size = 14)  +
  ylim(NA, 37) +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )


ggplot(pafip, aes(x = p3tto0, y = IncBMI3M)) +
  geom_boxplot(fill = "#69b3a2", outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(
    x = "Treatment",
    y = "BMI change in the first 3 months"
  ) +
  theme(
    axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  ) +
  ylim(NA, 5)

# ------------------------------
# 2. CLINICAL Linear Regression Models
# ------------------------------

# Variables
clinical_vars <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                   "sexo", "Edad", "p2diag6m", "Alcohol", "Tobacco", "Cannabis")

prs_vars <- c("BMI_prs", "T2D_prs", "BP_prs", "HYPCH_prs", "MET_SYN_prs", "NON_ISC_prs") #all PRSs

# Clinical model
clinical_model_BMI3M  <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars)))
clinical_model_BMI1a  <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars)))
clinical_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars))) 
clinical_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars)))
clinical_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars)))
clinical_model_BMI3a  <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars)))

summary(clinical_model_BMI3M)
summary(clinical_model_BMI1a)
summary(clinical_model_IncBMI3M)
summary(clinical_model_IncBMI1A)
summary(clinical_model_IncBMI3A)
summary(clinical_model_BMI3a)


# ------------------------------
# 3. CLINICAL + GENETIC Linear Regression Models (PRS Models)
# ------------------------------

# PRS BMI Model
prs_BMI_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BMI_prs))
prs_BMI_model_BMI3m <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), BMI_prs))
prs_BMI_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), BMI_prs))
prs_BMI_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), BMI_prs))
prs_BMI_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), BMI_prs))
prs_BMI_model_0 <- lm(BMI0 ~ ., data = pafip %>% select(BMI0, all_of(clinical_vars), BMI_prs))
prs_BMI_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), BMI_prs))

summary(prs_BMI_model_BMI3m)
summary(prs_BMI_model_BMI1a)
summary(prs_BMI_model_IncBMI3M)
summary(prs_BMI_model_IncBMI1A)
summary(prs_BMI_model_IncBMI3A)
summary(prs_BMI_model_0)
summary(prs_BMI_model_BMI3a)


# ALL PRS
prs_TOTAL_model_BMI3m <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), all_of(prs_vars)))
prs_TOTAL_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), all_of(prs_vars)))
prs_TOTAL_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), all_of(prs_vars)))
prs_TOTAL_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), all_of(prs_vars)))
prs_TOTAL_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), all_of(prs_vars)))
prs_TOTAL_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), all_of(prs_vars)))

summary(prs_TOTAL_model_BMI3m)
summary(prs_TOTAL_model_BMI1a)
summary(prs_TOTAL_model_IncBMI3M)
summary(prs_TOTAL_model_IncBMI1A)
summary(prs_TOTAL_model_IncBMI3A)
summary(prs_TOTAL_model_BMI3a)

# ALL PRS Except 1
prs_vars_except_1 <- prs_vars[prs_vars != "MET_SYN_prs"]
prs_TOTAL_model_BMI3m_except_1 <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_BMI1a_except_1 <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI3M_except_1 <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI1A_except_1 <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI3A_except_1 <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_BMI3a_except_1 <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), all_of(prs_vars_except_1)))

summary(prs_TOTAL_model_BMI3a_except_1)
summary(prs_TOTAL_model_BMI3m_except_1)
summary(prs_TOTAL_model_BMI1a_except_1)
summary(prs_TOTAL_model_IncBMI3M_except_1)
summary(prs_TOTAL_model_IncBMI1A_except_1)
summary(prs_TOTAL_model_IncBMI3A_except_1)

## Individual PRSs Models
# TYPE 2 DIABETES
prs_T2D_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), T2D_prs))
prs_T2D_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), T2D_prs))
prs_T2D_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), T2D_prs))
prs_T2D_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), T2D_prs))
prs_T2D_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), T2D_prs))
prs_T2D_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), T2D_prs))

summary(prs_T2D_model_BMI3a)
summary(prs_T2D_model_BMI3M)
summary(prs_T2D_model_BMI1a)
summary(prs_T2D_model_IncBMI3M)
summary(prs_T2D_model_IncBMI1A)
summary(prs_T2D_model_IncBMI3A)


# BLOOD PRESSURE
prs_BP_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), BP_prs))
prs_BP_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BP_prs))
prs_BP_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), BP_prs))
prs_BP_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), BP_prs))
prs_BP_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), BP_prs))
prs_BP_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), BP_prs))

summary(prs_BP_model_BMI3a)
summary(prs_BP_model_BMI3M)
summary(prs_BP_model_BMI1a)
summary(prs_BP_model_IncBMI3M)
summary(prs_BP_model_IncBMI1A)
summary(prs_BP_model_IncBMI3A)

# HYPCHOLESTEROLEMIA
prs_HYPCH_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), HYPCH_prs))
prs_HYPCH_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), HYPCH_prs))
prs_HYPCH_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), HYPCH_prs))
prs_HYPCH_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), HYPCH_prs))
prs_HYPCH_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), HYPCH_prs))
prs_HYPCH_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), HYPCH_prs))

summary(prs_HYPCH_model_BMI3a)
summary(prs_HYPCH_model_BMI3M)
summary(prs_HYPCH_model_BMI1a)
summary(prs_HYPCH_model_IncBMI3M)
summary(prs_HYPCH_model_IncBMI1A)
summary(prs_HYPCH_model_IncBMI3A)

# METABOLIC SYNDROME
prs_MET_SYN_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), MET_SYN_prs))
prs_MET_SYN_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), MET_SYN_prs))
prs_MET_SYN_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), MET_SYN_prs))
prs_MET_SYN_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), MET_SYN_prs))
prs_MET_SYN_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), MET_SYN_prs))
prs_MET_SYN_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), MET_SYN_prs))

summary(prs_MET_SYN_model_BMI3a)
summary(prs_MET_SYN_model_BMI3M)
summary(prs_MET_SYN_model_BMI1a)
summary(prs_MET_SYN_model_IncBMI3M)
summary(prs_MET_SYN_model_IncBMI1A)
summary(prs_MET_SYN_model_IncBMI3A)

# NON_ISCHEMIC CARDIOMYOPATHY
prs_NON_ISC_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), NON_ISC_prs))
prs_NON_ISC_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), NON_ISC_prs))
prs_NON_ISC_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), NON_ISC_prs))
prs_NON_ISC_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), NON_ISC_prs))
prs_NON_ISC_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), NON_ISC_prs))
prs_NON_ISC_model_BMI3a <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), NON_ISC_prs))

summary(prs_NON_ISC_model_BMI3a)
summary(prs_NON_ISC_model_BMI3M)
summary(prs_NON_ISC_model_BMI1a)
summary(prs_NON_ISC_model_IncBMI3M)
summary(prs_NON_ISC_model_IncBMI1A)
summary(prs_NON_ISC_model_IncBMI3A)


# PRS BMI + MET_SYN
prs_BMI_MET_SYN_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BMI_prs, MET_SYN_prs))
prs_BMI_MET_SYN_model_IncBMI3M <- lm(IncBMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(IncBMI3M, BMI0, all_of(clinical_vars), BMI_prs, MET_SYN_prs))
prs_BMI_MET_SYN_model_IncBMI1A <- lm(IncBMI1A ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(IncBMI1A, BMI0, all_of(clinical_vars), BMI_prs, MET_SYN_prs))

summary(prs_BMI_MET_SYN_model_BMI1a)
summary(prs_BMI_MET_SYN_model_IncBMI3M)
summary(prs_BMI_MET_SYN_model_IncBMI1A)

# COMBINATIONS
prs_NON_ISC_MET_model_BMI3M <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), NON_ISC_prs, MET_SYN_prs, BMI_prs))
prs_NON_ISC_MET_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), NON_ISC_prs, MET_SYN_prs, BMI_prs))
prs_NON_ISC_MET_model_IncBMI3M <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), NON_ISC_prs, MET_SYN_prs, BMI_prs))
prs_NON_ISC_MET_model_IncBMI1A <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), NON_ISC_prs, MET_SYN_prs, BMI_prs))
prs_NON_ISC_MET_model_IncBMI3A <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), NON_ISC_prs, MET_SYN_prs, BMI_prs))

summary(prs_NON_ISC_MET_model_BMI3M)
summary(prs_NON_ISC_MET_model_BMI1a)
summary(prs_NON_ISC_MET_model_IncBMI3M)
summary(prs_NON_ISC_MET_model_IncBMI1A)
summary(prs_NON_ISC_MET_model_IncBMI3A)



prs_BMI_model_BMI3m_ESP <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), BMI_prs, NON_ISC_prs))
prs_BMI_model_BMI1A_ESP <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BMI_prs, NON_ISC_prs))
prs_BMI_comb_BMI3m_ESP <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), BMI_prs, T2D_prs, MET_SYN_prs))

summary(prs_BMI_model_BMI3m_ESP)
summary(prs_BMI_model_BMI1A_ESP)
summary(prs_BMI_comb_BMI3m_ESP)


# ------------------------------
# 5. Compare models with ANOVA
# ------------------------------

# Clinical vs PRS_BMI
anova(clinical_model_BMI3M, prs_BMI_model_BMI3m)
anova(clinical_model_BMI1a, prs_BMI_model_BMI1a)
anova(clinical_model_IncBMI3M, prs_BMI_model_IncBMI3M)
anova(clinical_model_IncBMI1A, prs_BMI_model_IncBMI1A)
anova(clinical_model_IncBMI3A, prs_BMI_model_IncBMI3A)
anova(clinical_model_BMI3a, prs_BMI_model_BMI3a)

# PRS_BMI vs ALL_PRS
anova(prs_BMI_model_BMI3m, prs_TOTAL_model_BMI3m)
anova(prs_BMI_model_BMI1a, prs_TOTAL_model_BMI1a)
anova(prs_BMI_model_IncBMI3M, prs_TOTAL_model_IncBMI3M)
anova(prs_BMI_model_IncBMI1A, prs_TOTAL_model_IncBMI1A)
anova(prs_BMI_model_IncBMI3A, prs_TOTAL_model_IncBMI3A)
anova(prs_BMI_model_BMI3a, prs_TOTAL_model_BMI3a)

# Clinical vs Individual PRS
anova(clinical_model_BMI3M, prs_T2D_model_BMI3M)
anova(clinical_model_BMI3M, prs_BP_model_BMI3M)
anova(clinical_model_BMI3M, prs_HYPCH_model_BMI3M)
anova(clinical_model_BMI3M, prs_MET_SYN_model_BMI3M)
anova(clinical_model_BMI3M, prs_NON_ISC_model_BMI3M)

anova(clinical_model_BMI1a, prs_T2D_model_BMI1a)
anova(clinical_model_BMI1a, prs_BP_model_BMI1a)
anova(clinical_model_BMI1a, prs_HYPCH_model_BMI1a)
anova(clinical_model_BMI1a, prs_MET_SYN_model_BMI1a)
anova(clinical_model_BMI1a, prs_NON_ISC_model_BMI1a)

anova(clinical_model_IncBMI3M, prs_T2D_model_IncBMI3M)
anova(clinical_model_IncBMI3M, prs_BP_model_IncBMI3M)
anova(clinical_model_IncBMI3M, prs_HYPCH_model_IncBMI3M)
anova(clinical_model_IncBMI3M, prs_MET_SYN_model_IncBMI3M)
anova(clinical_model_IncBMI3M, prs_NON_ISC_model_IncBMI3M)

anova(clinical_model_IncBMI1A, prs_T2D_model_IncBMI1A)
anova(clinical_model_IncBMI1A, prs_BP_model_IncBMI1A)
anova(clinical_model_IncBMI1A, prs_HYPCH_model_IncBMI1A)
anova(clinical_model_IncBMI1A, prs_MET_SYN_model_IncBMI1A)
anova(clinical_model_IncBMI1A, prs_NON_ISC_model_IncBMI1A)

anova(clinical_model_IncBMI3A, prs_T2D_model_IncBMI3A)
anova(clinical_model_IncBMI3A, prs_BP_model_IncBMI3A)
anova(clinical_model_IncBMI3A, prs_HYPCH_model_IncBMI3A)
anova(clinical_model_IncBMI3A, prs_MET_SYN_model_IncBMI3A)
anova(clinical_model_IncBMI3A, prs_NON_ISC_model_IncBMI3A)

anova(clinical_model_BMI3a, prs_T2D_model_BMI3a)
anova(clinical_model_BMI3a, prs_BP_model_BMI3a)
anova(clinical_model_BMI3a, prs_HYPCH_model_BMI3a)
anova(clinical_model_BMI3a, prs_MET_SYN_model_BMI3a)
anova(clinical_model_BMI3a, prs_NON_ISC_model_BMI3a)

summary(clinical_model_IncBMI3A)
summary(prs_MET_SYN_model_IncBMI3A)

# Clinical vs All except BMI
anova(clinical_model_BMI3M, prs_TOTAL_model_BMI3m_except_1)
anova(clinical_model_BMI1a, prs_TOTAL_model_BMI1a_except_1)
anova(clinical_model_IncBMI3M, prs_TOTAL_model_IncBMI3M_except_1)
anova(clinical_model_IncBMI1A, prs_TOTAL_model_IncBMI1A_except_1)
anova(clinical_model_IncBMI3A, prs_TOTAL_model_IncBMI3A_except_1)

# PRS BMI vs All Except one

prs_vars_except_1 <- prs_vars[prs_vars != "BMI_prs"]
prs_TOTAL_model_BMI3m_except_1 <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_BMI1a_except_1 <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI3M_except_1 <- lm(INC3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(INC3M, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI1A_except_1 <- lm(INC1Y ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(INC1Y, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_IncBMI3A_except_1 <- lm(INC3Y ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(INC3Y, BMI0, all_of(clinical_vars), all_of(prs_vars_except_1)))
prs_TOTAL_model_BMI3a_except_1 <- lm(BMI3a ~ pafip$CPZ_1_3m + pafip$p3tto12 + ., data = pafip %>% select(BMI3a, all_of(clinical_vars), all_of(prs_vars_except_1)))


anova(prs_BMI_model_BMI3m, prs_TOTAL_model_BMI3m_except_1)
anova(prs_BMI_model_BMI1a, prs_TOTAL_model_BMI1a_except_1)
anova(prs_BMI_model_IncBMI3M, prs_TOTAL_model_IncBMI3M_except_1)
anova(prs_BMI_model_IncBMI1A, prs_TOTAL_model_IncBMI1A_except_1)
anova(prs_BMI_model_IncBMI3A, prs_TOTAL_model_IncBMI3A_except_1)

# Clinical vs All
anova(clinical_model_BMI3M, prs_TOTAL_model_BMI3m)
anova(clinical_model_BMI1a, prs_TOTAL_model_BMI1a)
anova(clinical_model_IncBMI3M, prs_TOTAL_model_IncBMI3M)
anova(clinical_model_IncBMI1A, prs_TOTAL_model_IncBMI1A)
anova(clinical_model_IncBMI3A, prs_TOTAL_model_IncBMI3A)
anova(clinical_model_BMI3a, prs_TOTAL_model_BMI3a)


# PRS BMI vs PRS BMI + MET_SYN

prs_BMI_METSYN_model_BMI3m <- lm(BMI3M ~ pafip$CPZ_0 + pafip$p3tto0 + ., data = pafip %>% select(BMI3M, all_of(clinical_vars), BMI_prs, MET_SYN_prs))
anova(prs_BMI_model_BMI3m, prs_BMI_METSYN_model_BMI3m)

prs_BMI_METSYN_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BMI_prs, MET_SYN_prs, NON_ISC_prs, T2D_prs))
anova(prs_BMI_model_BMI1a, prs_BMI_METSYN_model_BMI1a)
summary(prs_MET_SYN_model_BMI1a)

prs_BMI_NONISC_model_BMI1a <- lm(BMI1a ~ pafip$CPZ_1_3m + pafip$p3tto3m + ., data = pafip %>% select(BMI1a, all_of(clinical_vars), BMI_prs, NON_ISC_prs))
anova(prs_BMI_model_BMI1a, prs_BMI_NONISC_model_BMI1a)
anova(clinical_model_BMI1a, prs_BMI_NONISC_model_BMI1a)


summary(prs_BMI_NONISC_model_BMI1a)
summary(prs_BMI_model_BMI1a)
