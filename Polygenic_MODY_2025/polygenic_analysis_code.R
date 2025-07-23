####
# Polygenic MODY Analyisis
# Code used in Data Anlysis of "Common Genetic Variants Modify Disease Risk and Clinical Presentation in Monogenic Diabetes"
# Aim: Identify if polygenic risk influences presetnaion of HNF-MODY
####


rm(list=ls()) ##clear enviroment



#### load required libraries  ####
library(ggplot2)
library(ggpubr)
library(haven)
library(plyr)
library(tidyverse)

##for mixed effects
library(broom.mixed)
library(lme4)
library(performance)
library(lmerTest)





#### load stata datframe with phenotype and gentic risk score info ####
## this data frame has been cleaned, with predefined groups for T2D, Controls etc.
df <- read_dta('C:/Users/jm1493/OneDrive - University of Exeter/work/grs/local_topmed/analysis/stata_datasets/jml_mody_pheno.dta')


##define severe diabetes
df <- df %>%
  mutate(sevdiab = case_when(
    hbdcctkp >= 8.5 | ins_any == 1 ~ 1,
    hbdcctkp < 8.5 & ins_any == 0 ~ 0,
    TRUE ~ NA_real_                          
  ))



##main pgs scores
mscores <- c("stdt2dsuzuki24","stdt1d2", "prs_air_30", "prs_bmilocki",
             "prs_fg", "prs_hba1c", "prs_fi_80",
                 "stdlipodyslottaall_tphg38", "prs_whradjbmicombinedsexukbb")

##partitioned scores
pscores <- c("stdt2d_betacell_pineg","stdt2d_betacell_pipos", "stdt2d_bodyfat", "stdt2d_lipo" , "std_t2d_liverlipid",
                "stdt2d_metabollicsyn","stdt2d_obesity","stdt2d_residualglycaemic")


##all scores
scores <- c(mscores, pscores)


##labels for plotting
custom_labels <- c(
  "stdt2dsuzuki24" = "T2D",
  "stdt1d2" = "T1D",
  "prs_air_30" = "Acute Insulin Response",
  "prs_bmilocki" = "Body Mass Index",
  "prs_fg" = "Fasting Glucose",
  "prs_hba1c" = "HbA1c",
  "prs_fi_80" = "Fasting Insulin",
  "prs_whradjbmicombinedsexukbb" = "Waist Hip Ratio",
  "stdlipodyslottaall_tphg38" = " Lipodystrophy",
  "stdt2d_betacell_pineg" = "Beta Cell Proinsulin Neg",
  "stdt2d_betacell_pipos" = "Beta Cell Proinsulin Pos",
  "stdt2d_bodyfat" = "Body Fat",
  "stdt2d_lipo" = "Lipodystrophy",
  "std_t2d_liverlipid" = "Liver Lipid",
  "stdt2d_metabollicsyn" = "Metabollic Syndrome",
  "stdt2d_obesity" = "Obesity",
  "stdt2d_residualglycaemic" = "Residual Glycaemic")




####  Result 1: Raw Differnce from Controls ####

## unadjusted raw mean difference (standardised differnce)

##t2d and mody vs controls vs unknown
## and then per gene

##all mody carriers v controls

###loop through predefined groups

comparison_labels_1 <- c(
  "controlmodyall" ="MODY_CONTROL",
  "unkcontrol" ="UNSOLVED_CONTROL",
  "t2dcontrol" ="T2D_CONTROL",
  "modyvt2d" = "MODY_T2D",
  "modyvunknown" = "MODY_UNSOLVED",
  "unkvt2d" = "UNSOLVED_T2D",
  "controlcasebetacell" ="MODY_CONTROL_PROBANDS",
)


##varaibles to loop through
## e.g controls vs mody 
comparison_vars <- list("controlmodyall","t2dcontrol","unkcontrol", "modyvt2d", 
                        "modyvunknown", "unkvt2d", "controlcasebetacell")

results_df <- data.frame()

z_score <- 1.96

##scale pgs so controls are mean = 0, sd = 1

# Scaling controls to mean = 0, sd = 1
for (gprs in scores) {
  
  # Calculate mean and SD for the control group
  control_mean <- mean(df[[gprs]][df$controlmodyall == 0], na.rm = TRUE)
  control_sd <- sd(df[[gprs]][df$controlmodyall == 0], na.rm = TRUE)
  
  # Scale GRS for all individuals based on control mean and SD
  df[[paste0(gprs, "_scaled")]] <- (df[[gprs]] - control_mean) / control_sd
  

  ##now test the difference from controls for each pgs for each variable

  for (comp_var in comparison_vars) {
    
    # Create the formula for the linear model using the scaled GRS
    formula_str <- paste(paste0(gprs, "_scaled"), "~", comp_var, "+ eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10")
    formula <- as.formula(formula_str)
    
    # Fit the linear model
    model <- lm(formula, data = df)
    model_summary <- summary(model)
    coeff <- model_summary$coefficients[2, ] 
    
    estimate <- coeff[1]  
    se <- coeff[2]  
    p_value <- coeff[4]  
    
    # Calculate the GRS SD difference
    grs_sd_difference <- estimate  
    
    lower_ci <- grs_sd_difference - (z_score * se)
    upper_ci <- grs_sd_difference + (z_score * se)
    
    grs_label <- custom_labels[[gprs]]
    comparison_label <- comparison_labels_1[[comp_var]]
    
    results_df <- rbind(results_df, data.frame(
      PGS = grs_label,
      Comparison = comparison_label,
      PGS_SD_Difference = grs_sd_difference,
      SE = se,
      Lower_CI = lower_ci,
      Upper_CI = upper_ci,
      P_Value = p_value
    ))
  }
}


##now the the standradised differnce but per gene

results_df_gene <- data.frame()

comparison_vars <- list("controlmodyall")

valid_genes <- c("HNF1A", "HNF1B", "HNF4A")

for (gene_val in valid_genes) {
  
  
  filtered_data <- df[df$controlcasebetacell == 0 | (df$controlmodyall == 1 & df$gene == gene_val & df$diab01 == 1), ]
  
  for (gprs in scores) {
    
    
    control_mean <- mean(filtered_data[[gprs]][filtered_data$controlmodyall == 0], na.rm = TRUE)
    control_sd <- sd(filtered_data[[gprs]][filtered_data$controlmodyall == 0], na.rm = TRUE)
    
  
    filtered_data[[paste0(gprs, "_scaled")]] <- (filtered_data[[gprs]] - control_mean) / control_sd
    
  
    for (comp_var in comparison_vars) {
    
    
    
      
      
      formula_str <- paste(paste0(gprs, "_scaled"), "~", comp_var, "+ eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10")
      formula <- as.formula(formula_str)
      
      
      model <- lm(formula, data = filtered_data)
      model_summary <- summary(model)
      coeff <- model_summary$coefficients[2, ] 
      
      grs_sd_difference <- coeff[1]  
      se <- coeff[2]  
      p_value <- coeff[4]  
      
      
      lower_ci <- grs_sd_difference - (z_score * se)
      upper_ci <- grs_sd_difference + (z_score * se)
      
      
      grs_label <- custom_labels[[gprs]]
      
      
      results_df_gene <- rbind(results_df_gene, data.frame(
        PGS = grs_label,
        Comparison = gene_val,
        PGS_SD_Difference = grs_sd_difference,
        SE = se,
        Lower_CI = lower_ci,
        Upper_CI = upper_ci,
        P_Value = p_value
        
      ))
    }
  }
}


results_df_rawdiff <- rbind(results_df, results_df_gene)
  


## we next use logistic regression to determine which scores are independent
## mutliple formulae are used to detemine if indepent of other clincal characteristics

logistic_results_df <- data.frame()

# List of comparison variables and their labels
comparison_vars <- c("controlmodyall", "t2dcontrol", "unkcontrol")
comparison_labels <- c("MODY_CONTROL", "T2D_CONTROL", "UNSOLVED_CONTROL")

# List of adjustment formulas and their corresponding labels
full_adjustment <- "stdt2dsuzuki24_scaled + stdt1d2_scaled + prs_bmilocki_scaled + prs_fg_scaled + prs_hba1c_scaled + prs_air_30_scaled + prs_whradjbmicombinedsexukbb_scaled + prs_fi_80_scaled + stdlipodyslottaall_tphg38_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10"
full_adjustment_clusters <- "stdt2d_betacell_pineg_scaled + stdt2d_betacell_pipos_scaled + stdt2d_bodyfat_scaled + stdt2d_lipo_scaled + stdt2d_metabollicsyn_scaled + stdt2d_obesity_scaled + stdt2d_residualglycaemic_scaled + std_t2d_liverlipid_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10"

full_adjustment_withclin <- "stdt2dsuzuki24_scaled + stdt1d2_scaled + prs_bmilocki_scaled + prs_fg_scaled + prs_hba1c_scaled + prs_air_30_scaled + prs_whradjbmicombinedsexukbb_scaled + prs_fi_80_scaled + stdlipodyslottaall_tphg38_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +sexall +bmiadultextdare + age"
full_adjustment_clusters_withclin <- "stdt2d_betacell_pineg_scaled + stdt2d_betacell_pipos_scaled + stdt2d_bodyfat_scaled + stdt2d_lipo_scaled + stdt2d_metabollicsyn_scaled + stdt2d_obesity_scaled + stdt2d_residualglycaemic_scaled + std_t2d_liverlipid_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +sexall +bmiadultextdare + age"

full_adjustment_withclinmf <- "stdt2dsuzuki24_scaled + stdt1d2_scaled + prs_bmilocki_scaled + prs_fg_scaled + prs_hba1c_scaled + prs_air_30_scaled + prs_whradjbmicombinedsexukbb_scaled + prs_fi_80_scaled + stdlipodyslottaall_tphg38_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +sexall +bmiadultextdare + age +mf"
full_adjustment_clusters_withclinmf <- "stdt2d_betacell_pineg_scaled + stdt2d_betacell_pipos_scaled + stdt2d_bodyfat_scaled + stdt2d_lipo_scaled + stdt2d_metabollicsyn_scaled + stdt2d_obesity_scaled + stdt2d_residualglycaemic_scaled + std_t2d_liverlipid_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +sexall +bmiadultextdare + age +mf"

full_adjustment_mf <- "stdt2dsuzuki24_scaled + stdt1d2_scaled + prs_bmilocki_scaled + prs_fg_scaled + prs_hba1c_scaled + prs_air_30_scaled + prs_whradjbmicombinedsexukbb_scaled + prs_fi_80_scaled + stdlipodyslottaall_tphg38_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +mf"
full_adjustment_clusters_mf <- "stdt2d_betacell_pineg_scaled + stdt2d_betacell_pipos_scaled + stdt2d_bodyfat_scaled + stdt2d_lipo_scaled + stdt2d_metabollicsyn_scaled + stdt2d_obesity_scaled + stdt2d_residualglycaemic_scaled + std_t2d_liverlipid_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +mf"

full_adjustment_withclin_bw <- "stdt2dsuzuki24_scaled + stdt1d2_scaled + prs_bmilocki_scaled + prs_fg_scaled + prs_hba1c_scaled + prs_air_30_scaled + prs_whradjbmicombinedsexukbb_scaled + prs_fi_80_scaled + stdlipodyslottaall_tphg38_scaled + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +sexall +bmiadultextdare + age + birthweightg"

adjustment_formulas <- list(
  full_adjustment,
  full_adjustment_clusters,
  full_adjustment_withclin,
  full_adjustment_clusters_withclin,
  full_adjustment_withclinmf,
  full_adjustment_clusters_withclinmf,
  full_adjustment_mf,
  full_adjustment_clusters_mf,
  full_adjustment_withclin_bw
)

adjustment_labels <- c(
  "Full Adjustment",
  "Full Adjustment Clusters",
  "Full Adjustment with Clinical",
  "Full Adjustment Clusters with Clinical",
  "Full Adjustment with Clinical + MF",
  "Full Adjustment Clusters with Clinical + MF",
  "Full Adjustment + MF",
  "Full Adjustment Clusters + MF",
  "Full Adjustment with Clinical + Birth Weight"
)

# Loop through each comparison variable
for (comp_var in comparison_vars) {
  
  
  for (j in seq_along(adjustment_formulas)) {
    adjustment_formula <- adjustment_formulas[[j]]
    adjustment_label <- adjustment_labels[j]
    
    # Construct formula
    formula_str <- paste(comp_var, "~", adjustment_formula)
    formula <- as.formula(formula_str)
    
    # Fit logistic model
    logistic_model <- glm(formula, data = df, family = binomial)
    model_summary <- summary(logistic_model)
    
    
    for (i in 1:nrow(model_summary$coefficients)) {
      coeff <- model_summary$coefficients[i, ]
      
      log_odds_estimate <- coeff[1]
      se <- coeff[2]
      p_value <- coeff[4]
      
      variable_name <- rownames(model_summary$coefficients)[i]
      
      if (variable_name != "(Intercept)") {
        standardized_estimate <- log_odds_estimate
        
        # Odds ratio calculation
        odds_ratio_estimate <- exp(standardized_estimate)
        
        # Confidence intervals for log odds and odds ratio
        ci_lower_log_odds <- standardized_estimate - (z_score * se)
        ci_upper_log_odds <- standardized_estimate + (z_score * se)
        ci_lower_odds_ratio <- exp(ci_lower_log_odds)
        ci_upper_odds_ratio <- exp(ci_upper_log_odds)
        
     
        logistic_results_df <- rbind(logistic_results_df, data.frame(
          PGS = variable_name,
          Comparison = comparison_labels[which(comparison_vars == comp_var)],
          Adjustment_Label = adjustment_label,  
          Estimate_Log_Odds = standardized_estimate,
          Estimate_Odds_Ratio = odds_ratio_estimate,
          CI_Lower = ci_lower_odds_ratio,
          CI_Upper = ci_upper_odds_ratio,
          SE = se,
          P_Value = p_value
        ))
      }
    }
  }
}

logistic_results_df$PGS <- gsub("_scaled", "", logistic_results_df$PGS)
logistic_results_df <- logistic_results_df %>% filter(PGS %in% scores)
logistic_results_df <- read.xlsx("independent_scores.xlsx")
#write.xlsx(logistic_results_df, "independent_scores.xlsx")  



## t2d risk score before and after adjusting family history

results_df_mf <- data.frame(
  PGS = character(),
  Comparison = character(),
  Adjustment = character(),
  PGS_SD_Difference = numeric(),
  SE = numeric(),
  Lower_CI = numeric(),
  Upper_CI = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

comparison_vars <- list("controlmodyall","t2dcontrol","unkcontrol")
nscores <- list("stdt2dsuzuki24_scaled")

for (comp_var in comparison_vars) {
  for (gprs in nscores) {
    
    # Model without adjusting for 'mf'
    formula_str_no_mf <- paste(gprs, "~", comp_var, "+ eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10")
    formula_no_mf <- as.formula(formula_str_no_mf)
    
    model_no_mf <- lm(formula_no_mf, data = df)
    model_summary_no_mf <- summary(model_no_mf)
    coeff_no_mf <- model_summary_no_mf$coefficients[2, ]  # Get the coefficient for the comparison variable
    
    grs_sd_difference_no_mf <- coeff_no_mf[1]  
    se_no_mf <- coeff_no_mf[2]
    p_value_no_mf <- coeff_no_mf[4] 
  
    
    lower_ci_no_mf <- grs_sd_difference_no_mf - (z_score * se_no_mf)
    upper_ci_no_mf <- grs_sd_difference_no_mf + (z_score * se_no_mf)
    
    formula_str_with_mf <- paste(gprs, "~", comp_var, "+ eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + mf")
    formula_with_mf <- as.formula(formula_str_with_mf)
    
    model_with_mf <- lm(formula_with_mf, data = df)
    model_summary_with_mf <- summary(model_with_mf)
    coeff_with_mf <- model_summary_with_mf$coefficients[2, ]  # Get the coefficient for the comparison variable
    
    grs_sd_difference_with_mf <- coeff_with_mf[1]  # Estimate
    se_with_mf <- coeff_with_mf[2]  # Standard error
    p_value_with_mf <- coeff_with_mf[4]  # P-value
    
    
    lower_ci_with_mf <- grs_sd_difference_with_mf - (z_score * se_with_mf)
    upper_ci_with_mf <- grs_sd_difference_with_mf + (z_score * se_with_mf)
    
    # Labels for the current PGS and comparison variable
    grs_label <- "T2D"
    comparison_label <- comparison_labels_1[[comp_var]]
    
    # Add results before adjusting for 'mf'
    results_df_mf <- rbind(results_df_mf, data.frame(
      PGS = grs_label,
      Comparison = comparison_label,
      Adjustment = "Unadjusted",
      PGS_SD_Difference = grs_sd_difference_no_mf,
      SE = se_no_mf,
      Lower_CI = lower_ci_no_mf,
      Upper_CI = upper_ci_no_mf,
      P_Value = p_value_no_mf
    ))
    
    # Add results after adjusting for 'mf'
    results_df_mf <- rbind(results_df_mf, data.frame(
      PGS = grs_label,
      Comparison = comparison_label,
      Adjustment = "Adjusted for mf",
      PGS_SD_Difference = grs_sd_difference_with_mf,
      SE = se_with_mf,
      Lower_CI = lower_ci_with_mf,
      Upper_CI = upper_ci_with_mf,
      P_Value = p_value_with_mf
    ))
  }
}






## adjuseted mixed effect models --> effect on age diag, diabetes severity   ####

##  age of diagnosis



regagediag_df <- data.frame()
comparison_vars <- list("controlmodyall", "t2dcontrol", "unkcontrol", "controlcasebetacell", "hnf1a_control",
                        "hnf4a_control","hnf1b_control")
outcome_var <- "agediag"

for (score_var in scores) {
  
  ## standardise score so controls mean 0 and sd 1
  

  
  df <- df %>%
    mutate(standardized_score = (.data[[score_var]]))
  
  for (regroup_val in comparison_vars) {
    
    # Filter data based on category
    filtered_data <- df[!is.na(df[[regroup_val]]), ]
    
    
    
    if (regroup_val == "controlmodyall") {
      filtered_data <- filtered_data %>% filter(controlmodyall == 1)
      group_name <- "HNF_MODY"
    } else if (regroup_val == "t2dcontrol") {
      filtered_data <- filtered_data %>% filter(t2dcontrol == 1)
      group_name <- "T2D"
    } else if (regroup_val == "unkcontrol") {
      filtered_data <- filtered_data %>% filter(unkcontrol == 1)
      group_name <- "Unsolved" }
      else if (regroup_val %in% c("hnf1a_control","hnf4a_control","hnf1b_control")) {
        filtered_data <- filtered_data %>% filter(.data[[regroup_val]] == 1)
        group_name <- regroup_val
    } else if (regroup_val == "controlcasebetacell") {
      filtered_data <- filtered_data %>% filter(controlcasebetacell == 1)
      group_name <- "HNF_MODY_Probands"
    }
    

    
    # Define model formulas
    # Use Adjusted and Unadjusted models to see if effects independent

    if (regroup_val %in% c("controlmodyall")) {
      # Use mixed effects models for family-based data
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score +  eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + (1 | familyid2)"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + modyclin + exonintronaff + gene0 + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + bmiadultextdare + birthweightg + (1 | familyid2)"))
    } else if (regroup_val %in% c("hnf1a_control","hnf4a_control","hnf1b_control")) {
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score +  eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + (1 | familyid2)"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + modyclin + exonintronaff + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + bmiadultextdare  + (1 | familyid2)"))
    }else if (regroup_val %in% c("controlcasebetacell")) {
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + exonintronaff + gene0 + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + bmiadultextdare"))
    } else {
 
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + bmiadultextdare"))
    }
    
    # Fit unadjusted model and extract 
    if (regroup_val %in% c("controlmodyall", "hnf1a_control","hnf4a_control","hnf1b_control")) {
      model_unadj <- lmer(unadj_formula, data = filtered_data, REML = FALSE)
      model_summary_unadj <- broom.mixed::tidy(model_unadj, effects = "fixed", conf.int = TRUE)
    } else {
      model_unadj <- lm(unadj_formula, data = filtered_data)
      model_summary_unadj <- broom::tidy(model_unadj, conf.int = TRUE)
    }
    
    unadj_row <- model_summary_unadj %>%
      filter(term == "standardized_score") %>%
      summarise(
        outcome = outcome_var,
        term = score_var,
        comparison = regroup_val,
        fixed_effects = estimate,
        std_error = std.error,
        lower_ci = conf.low,
        upper_ci = conf.high,
        p_value = p.value,
        adj = "no"
      )
    
    # Fit adjusted model and extract 
    if (regroup_val %in% c("controlmodyall", "hnf1a_control","hnf4a_control","hnf1b_control")) {
      model_adj <- lmer(adj_formula, data = filtered_data, REML = FALSE)
      model_summary_adj <- broom.mixed::tidy(model_adj, effects = "fixed", conf.int = TRUE)
    } else {
      model_adj <- lm(adj_formula, data = filtered_data)
      model_summary_adj <- broom::tidy(model_adj, conf.int = TRUE)
    }
    
    adj_row <- model_summary_adj %>%
      filter(term == "standardized_score") %>%
      summarise(
        outcome = outcome_var,
        term = score_var,
        comparison = regroup_val,
        fixed_effects = estimate,
        std_error = std.error,
        lower_ci = conf.low,
        upper_ci = conf.high,
        p_value = p.value,
        adj = "yes"
      )
    

    regagediag_df <- bind_rows(regagediag_df, unadj_row, adj_row)
  }
}




####severe diabetes loop####

severe_df <- data.frame()
comparison_vars <- list("controlmodyall", "t2dcontrol", "unkcontrol", "controlcasebetacell", "hnf1a_control",
                        "hnf4a_control","hnf1b_control")
outcome_var <- "sevdiab"

for (score_var in scores) {
  
  
  df <- df %>%
    mutate(standardized_score = (.data[[score_var]]))
  
  for (regroup_val in comparison_vars) {
    
    # Filter data based on category
    filtered_data <- df[!is.na(df[[regroup_val]]), ]
    
    if (regroup_val == "controlmodyall") {
      filtered_data <- filtered_data %>% filter(controlmodyall == 1)
      group_name <- "HNF_MODY"
    } else if (regroup_val == "t2dcontrol") {
      filtered_data <- filtered_data %>% filter(t2dcontrol == 1)
      group_name <- "T2D"
    } else if (regroup_val == "unkcontrol") {
      filtered_data <- filtered_data %>% filter(unkcontrol == 1)
      group_name <- "Unsolved"
     } else if (regroup_val %in% c("hnf1a_control","hnf4a_control","hnf1b_control")) {
        filtered_data <- filtered_data %>% filter(.data[[regroup_val]] == 1)
        group_name <- regroup_val
    } else if (regroup_val == "controlcasebetacell") {
      filtered_data <- filtered_data %>% filter(controlcasebetacell == 1)
      group_name <- "HNF_MODY_Probands"
    }
    
    # Define model formulas
    if (regroup_val %in% c("controlmodyall")) {
      # Use mixed effects models for family-based data (controlmodyall)
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 +(1 | familyid2)"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + modyclin + exonintronaff + gene0 + 
                                     eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + 
                                     bmiadultextdare + bweightg + (1 | familyid2)"))
    } else if (regroup_val %in% c("hnf1a_control","hnf4a_control","hnf1b_control")) {
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score +  eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + (1 | familyid2)"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + modyclin + exonintronaff  + 
                                     eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + 
                                     bmiadultextdare + bweightg + (1 | familyid2)"))
    }
    
    else if (regroup_val %in% c("controlcasebetacell")) {
      # Use mixed effects models for family-based data (controlcasebetacell)
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score "))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + mf + yeardiag + exonintronaff + gene0 + 
                                     eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + 
                                     bmiadultextdare  "))
    } else {
      # Use standard models for non-family data
      unadj_formula <- as.formula(paste(outcome_var, "~ standardized_score"))
      adj_formula <- as.formula(paste(outcome_var, "~ standardized_score + sexall + eupc1 + eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + bmiadultextdare"))
    }
    
    # Fit unadjusted model and extract 
    if (regroup_val %in% c("controlmodyall","hnf1a_control","hnf4a_control","hnf1b_control")) {
      model_unadj <- glmer(unadj_formula, data = filtered_data, family = binomial(link = "logit"))
      model_summary_unadj <- broom.mixed::tidy(model_unadj, effects = "fixed", conf.int = TRUE)
    } else {
      model_unadj <- glm(unadj_formula, data = filtered_data, family = binomial(link = "logit"))
      model_summary_unadj <- broom::tidy(model_unadj, conf.int = TRUE)
    }
    
    unadj_row <- model_summary_unadj %>%
      filter(term == "standardized_score") %>%
      summarise(
        outcome = outcome_var,
        term = score_var,
        comparison = regroup_val,
        fixed_effects = estimate,
        std_error = std.error,
        lower_ci = conf.low,
        upper_ci = conf.high,
        p_value = p.value,
        adj = "no"
      )
    
    # Fit adjusted model and extract 
    if (regroup_val %in% c("controlmodyall","hnf1a_control","hnf4a_control","hnf1b_control")) {
      model_adj <- glmer(adj_formula, data = filtered_data, family = binomial(link = "logit"))
      model_summary_adj <- broom.mixed::tidy(model_adj, effects = "fixed", conf.int = TRUE)
    } else {
      model_adj <- glm(adj_formula, data = filtered_data, family = binomial(link = "logit"))
      model_summary_adj <- broom::tidy(model_adj, conf.int = TRUE)
    }
    
    adj_row <- model_summary_adj %>%
      filter(term == "standardized_score") %>%
      summarise(
        outcome = outcome_var,
        term = score_var,
        comparison = regroup_val,
        fixed_effects = estimate,
        std_error = std.error,
        lower_ci = conf.low,
        upper_ci = conf.high,
        p_value = p.value,
        adj = "yes"
      )
    
    # Combine unadjusted and adjusted results
    severe_df <- bind_rows(severe_df, unadj_row, adj_row)
  }
}

severe_df$odds_ratio <- exp(severe_df$fixed_effects) 
severe_df$odds_ratio_lower_ci <- exp(severe_df$lower_ci)
severe_df$odds_ratio_upper_ci <- exp(severe_df$upper_ci)
severe_df$or_se = exp(severe_df$fixed_effects) * severe_df$std_error

##next combine all risk scores into one model to see which are indpeendent
## mixed model adjusting score for other scores to see which independent 

adj_formula <- as.formula(paste("agediag", "~ stdt2dsuzuki24_standardized + stdt1d2_standardized + prs_air_30_standardized + prs_bmilocki_standardized +
              prs_fg_standardized + prs_hba1c_standardized + prs_fi_80_standardized +
              stdlipodyslottaall_tphg38_standardized + prs_whradjbmicombinedsexukbb_standardized + eupc1 + 
              eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + (1| familyid)"))


filtered_data <- df %>% filter(controlmodyall ==1)
model_adj <- lmer(adj_formula, data = filtered_data, REML = FALSE)

# Get the summary of the model
model_summary_adj <- broom.mixed::tidy(model_adj, effects = "fixed", conf.int = TRUE)
model_summary_adj$term <- gsub("_standardized$", "", model_summary_adj$term)
view(model_summary_adj)


##  mixed model adjusting for other scores (Diabetes Severity)
adj_formula <- as.formula(paste("sevdiab", "~ stdt2dsuzuki24_standardized + stdt1d2_standardized + prs_air_30_standardized + prs_bmilocki_standardized +
              prs_fg_standardized + prs_hba1c_standardized + prs_fi_80_standardized +
              stdlipodyslottaall_tphg38_standardized + prs_whradjbmicombinedsexukbb_standardized + eupc1 + 
              eupc2 + eupc3 + eupc4 + eupc5 + eupc6 + eupc7 + eupc8 + eupc9 + eupc10 + (1| familyid)"))
model_adj_logistic <- glmer(adj_formula, data = filtered_data, 
                            family = binomial(link = "logit"))
model_summary_logistic <- broom.mixed::tidy(model_adj_logistic, effects = "fixed", conf.int = TRUE)
model_summary_logistic$odds_ratio <- exp(model_summary_logistic$estimate) 
model_summary_logistic$odds_ratio_lower_ci <- exp(model_summary_logistic$conf.low)
model_summary_logistic$odds_ratio_upper_ci <- exp(model_summary_logistic$conf.high)
model_summary_logistic$term <- gsub("_standardized$", "", model_summary_logistic$term)



##get parental history data
parental_df <-  df %>%
  filter(!is.na(mf) & !is.na(group_var)) %>% ##mf == parental history (0 =none, 1 = 1 parent, 2 = both)  group_var == T2D, MODY etc.
  group_by(group_var, mf) %>%
  summarize(count = n(), .groups = "drop") %>%
  group_by(group_var) %>%
  mutate(percent = (count / sum(count)) * 100) %>%
  mutate(
    ci_lower_par = sapply(count, function(x) prop.test(x, sum(count))$conf.int[1]) * 100,  # Multiply by 100 to convert to percentage
    ci_upper_par = sapply(count, function(x) prop.test(x, sum(count))$conf.int[2]) * 100   # Multiply by 100 to convert to percentage
  ) %>%
  mutate(group_var = factor(group_var, levels = c("Control", "T2D", "Unsolved","MODY"))) %>%
  mutate(mf = factor(mf, levels = c(2, 1, 0))) %>%
  select(group_var, mf, percent)




##Compare T2D risk across difference variant groupings

##missense versus ptv


##capture summary info
# Genes including 'All'

genes_of_interest <- c("All", "HNF1A", "HNF1B", "HNF4A")


summary_table <- data.frame()

# Loop through each gene (including 'All')
for (gene in genes_of_interest) {
  
  df_filtered <- df %>%
    filter(!is.na(mis_ptv), !is.na(stdt2dsuzuki24_scaled))  ##mis ptv, split by revel <0.9, revel >=0.9, and ptvs
  
  if (gene != "All") {
    df_filtered <- df_filtered %>% filter(gene == !!gene)
  }

  # Calculate group stats
  summary_df <- df_filtered %>%
    group_by(mis_ptv) %>%
    summarise(
      mean_t2d_risk = mean(stdt2dsuzuki24_scaled, na.rm = TRUE),
      sd = sd(stdt2dsuzuki24_scaled, na.rm = TRUE),
      n = n(),
      se = sd / sqrt(n),
      .groups = "drop"
    )

  # append ANOVA p values
  if (n_distinct(df_filtered$mis_ptv) > 1) {
    anova_result <- aov(stdt2dsuzuki24_scaled ~ mis_ptv, data = df_filtered)
    anova_p_value <- summary(anova_result)[[1]]["Pr(>F)"][[1]][1]
  } else {
    anova_p_value <- NA  
  }

  # Add gene name and anova p to each row
  summary_df <- summary_df %>%
    mutate(
      gene = gene,
      anova_p_value = anova_p_value
    ) %>%
    select(gene, mis_ptv, mean_t2d_risk, se, anova_p_value)

  summary_table <- bind_rows(summary_table, summary_df)
}



