packages <- c(
  "survival",
  "survminer",
  "dplyr",
  "ggplot2",
  "stringr",
  "tidyr",
  "readr"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}


# ============================================================
# 1. Settings
# ============================================================

patient_file  <- "../metadata/data_clinical_patient.txt"
sample_file   <- "../metadata/data_clinical_sample.txt"
timeline_file <- "../metadata/data_timeline_status.txt"

output_dir <- "survival_results_final"
dir.create(output_dir, showWarnings = FALSE)

X_AXIS_MAX <- 35

SHOW_CONFIDENCE_INTERVALS <- FALSE

# If FALSE: p-value uses full data
ADMINISTRATIVE_CENSORING <- FALSE

# ============================================================
# 2. Read files
# ============================================================

read_file <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    comment.char = "#",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

patient <- read_file(patient_file)
sample <- read_file(sample_file)

timeline <- read.delim(
  timeline_file,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# ============================================================
# 3. Select patient-level survival metadata
# ============================================================

patient_surv <- patient %>%
  select(
    PATIENT_ID,
    OS_STATUS,
    OS_MONTHS,
    DSS_STATUS,
    DSS_MONTHS,
    PFS_STATUS,
    PFS_MONTHS,
    AGE,
    SEX,
    AJCC_PATHOLOGIC_TUMOR_STAGE
  )


# ============================================================
# 4. Select sample-level group metadata
# ============================================================

sample_group <- sample %>%
  select(
    PATIENT_ID,
    SAMPLE_ID,
    SAMPLE_TYPE
  ) %>%
  filter(SAMPLE_TYPE %in% c("Primary", "Metastasis"))


# ============================================================
# 5. Extract Submitted Specimen Dx date
# ============================================================
# START_DATE is in days from initial pathologic diagnosis.
# We convert it to months and use it as the sample-based time origin.
# If a patient has multiple Submitted Specimen Dx records, we take the earliest.

specimen_dx <- timeline %>%
  filter(STATUS == "Submitted Specimen Dx") %>%
  mutate(
    START_DATE = as.numeric(START_DATE),
    SAMPLE_DX_MONTHS = START_DATE / 30.44
  ) %>%
  filter(!is.na(SAMPLE_DX_MONTHS)) %>%
  group_by(PATIENT_ID) %>%
  summarise(
    SAMPLE_DX_DAYS = min(START_DATE, na.rm = TRUE),
    SAMPLE_DX_MONTHS = min(SAMPLE_DX_MONTHS, na.rm = TRUE),
    N_SPECIMEN_DX_RECORDS = n(),
    .groups = "drop"
  )


# ============================================================
# 6. Merge sample type + survival + sample diagnosis date
# ============================================================

df <- sample_group %>%
  inner_join(patient_surv, by = "PATIENT_ID") %>%
  left_join(specimen_dx, by = "PATIENT_ID")


# ============================================================
# 7. Prepare survival variables
# ============================================================
# OS_STATUS:  1:DECEASED / 0:LIVING
# DSS_STATUS: 1:DEAD WITH TUMOR / 0:ALIVE OR DEAD TUMOR FREE
# PFS_STATUS: 1:PROGRESSION / 0:CENSORED

df <- df %>%
  mutate(
    OS_MONTHS  = as.numeric(OS_MONTHS),
    DSS_MONTHS = as.numeric(DSS_MONTHS),
    PFS_MONTHS = as.numeric(PFS_MONTHS),
    
    OS_EVENT  = ifelse(str_detect(as.character(OS_STATUS), "^1:"), 1, 0),
    DSS_EVENT = ifelse(str_detect(as.character(DSS_STATUS), "^1:"), 1, 0),
    PFS_EVENT = ifelse(str_detect(as.character(PFS_STATUS), "^1:"), 1, 0),
    
    SAMPLE_TYPE = factor(SAMPLE_TYPE, levels = c("Metastasis", "Primary")),
    SEX = factor(SEX),
    AJCC_PATHOLOGIC_TUMOR_STAGE = factor(AJCC_PATHOLOGIC_TUMOR_STAGE)
  )


# ============================================================
# 8. Correct survival time from sample diagnosis
# ============================================================

df <- df %>%
  mutate(
    OS_MONTHS_FROM_SAMPLE  = OS_MONTHS  - SAMPLE_DX_MONTHS,
    DSS_MONTHS_FROM_SAMPLE = DSS_MONTHS - SAMPLE_DX_MONTHS,
    PFS_MONTHS_FROM_SAMPLE = PFS_MONTHS - SAMPLE_DX_MONTHS,
    
    OS_VALID_FROM_SAMPLE  = !is.na(OS_MONTHS_FROM_SAMPLE)  & OS_MONTHS_FROM_SAMPLE  > 0,
    DSS_VALID_FROM_SAMPLE = !is.na(DSS_MONTHS_FROM_SAMPLE) & DSS_MONTHS_FROM_SAMPLE > 0,
    PFS_VALID_FROM_SAMPLE = !is.na(PFS_MONTHS_FROM_SAMPLE) & PFS_MONTHS_FROM_SAMPLE > 0
  )


# ============================================================
# 9. Save merged dataset
# ============================================================

write.csv(
  df,
  file = file.path(output_dir, "merged_survival_primary_metastasis_time_corrected.csv"),
  row.names = FALSE
)


# ============================================================
# 10. Save data summary
# ============================================================

summary_file <- file.path(output_dir, "data_summary.txt")

sink(summary_file)

cat("============================================================\n")
cat("TCGA-SKCM survival analysis summary\n")
cat("============================================================\n\n")

cat("Original sample type counts:\n")
print(table(df$SAMPLE_TYPE, useNA = "ifany"))

cat("\nSubmitted Specimen Dx records per patient:\n")
print(summary(df$N_SPECIMEN_DX_RECORDS))

cat("\nSample diagnosis months summary:\n")
print(summary(df$SAMPLE_DX_MONTHS))

cat("\nOS months from initial diagnosis:\n")
print(summary(df$OS_MONTHS))

cat("\nOS months from sample diagnosis:\n")
print(summary(df$OS_MONTHS_FROM_SAMPLE))

cat("\nDSS months from initial diagnosis:\n")
print(summary(df$DSS_MONTHS))

cat("\nDSS months from sample diagnosis:\n")
print(summary(df$DSS_MONTHS_FROM_SAMPLE))

cat("\nPFS months from initial diagnosis:\n")
print(summary(df$PFS_MONTHS))

cat("\nPFS months from sample diagnosis:\n")
print(summary(df$PFS_MONTHS_FROM_SAMPLE))

cat("\nValid OS from sample counts:\n")
print(table(df$SAMPLE_TYPE[df$OS_VALID_FROM_SAMPLE], useNA = "ifany"))

cat("\nValid DSS from sample counts:\n")
print(table(df$SAMPLE_TYPE[df$DSS_VALID_FROM_SAMPLE], useNA = "ifany"))

cat("\nValid PFS from sample counts:\n")
print(table(df$SAMPLE_TYPE[df$PFS_VALID_FROM_SAMPLE], useNA = "ifany"))

cat("\nOS event counts by sample type:\n")
print(table(df$SAMPLE_TYPE, df$OS_EVENT, useNA = "ifany"))

cat("\nDSS event counts by sample type:\n")
print(table(df$SAMPLE_TYPE, df$DSS_EVENT, useNA = "ifany"))

cat("\nPFS event counts by sample type:\n")
print(table(df$SAMPLE_TYPE, df$PFS_EVENT, useNA = "ifany"))

sink()

cat("\nSummary saved to:", summary_file, "\n")


# ============================================================
# 11. Helper function for p-value label
# ============================================================

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("Logrank Test P-Value: NA")
  }
  
  if (p_value < 0.0001) {
    return("Logrank Test P-Value: < 0.0001")
  }
  
  if (p_value < 0.001) {
    return(paste0(
      "Logrank Test P-Value: ",
      formatC(p_value, format = "e", digits = 2)
    ))
  }
  
  return(paste0(
    "Logrank Test P-Value: ",
    signif(p_value, digits = 3)
  ))
}


# ============================================================
# 12. Kaplan-Meier plotting function
# ============================================================

plot_km_cbio_style <- function(data,
                               time_col,
                               event_col,
                               title,
                               x_label,
                               y_label,
                               output_name,
                               xlim_max = 60,
                               break_time_by = 5,
                               administrative_censoring = FALSE,
                               show_ci = FALSE) {
  
  km_df <- data %>%
    filter(
      !is.na(.data[[time_col]]),
      !is.na(.data[[event_col]]),
      !is.na(SAMPLE_TYPE),
      .data[[time_col]] > 0
    ) %>%
    mutate(
      original_time = .data[[time_col]],
      original_event = .data[[event_col]],
      
      KM_GROUP = factor(
        ifelse(SAMPLE_TYPE == "Metastasis", "(A) Metastasis", "(B) Primary"),
        levels = c("(A) Metastasis", "(B) Primary")
      )
    )
  
  if (administrative_censoring) {
    km_df <- km_df %>%
      mutate(
        surv_time = pmin(original_time, xlim_max),
        surv_event = ifelse(original_time <= xlim_max, original_event, 0)
      )
  } else {
    km_df <- km_df %>%
      mutate(
        surv_time = original_time,
        surv_event = original_event
      )
  }
  
  if (length(unique(km_df$KM_GROUP)) < 2) {
    cat("Skipping plot: fewer than two groups available.\n")
    return(NULL)
  }
  
  fit <- survfit(
    Surv(surv_time, surv_event) ~ KM_GROUP,
    data = km_df
  )
  
  logrank <- survdiff(
    Surv(surv_time, surv_event) ~ KM_GROUP,
    data = km_df
  )
  
  p_value <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  p_label <- format_p_value(p_value)
  
  p <- ggsurvplot(
    fit,
    data = km_df,
    
    conf.int = show_ci,
    risk.table = TRUE,
    risk.table.height = 0.22,
    
    censor = TRUE,
    censor.shape = "+",
    censor.size = 3.2,
    
    palette = c("#E64B35", "#2E86DE"),
    
    legend.title = "",
    legend.labs = c("(A) Metastasis", "(B) Primary"),
    legend = "right",
    
    xlim = c(0, xlim_max),
    break.time.by = break_time_by,
    
    surv.scale = "percent",
    
    xlab = paste0(x_label, " (Months)"),
    ylab = paste0("Probability of ", y_label),
    
    pval = FALSE,
    
    risk.table.y.text = TRUE,
    risk.table.y.text.col = FALSE,
    risk.table.title = "Number at risk (n)",
    risk.table.fontsize = 4.5,
    
    ggtheme = theme_classic(base_size = 14)
  )
  
  p$plot <- p$plot +
    ggtitle(title) +
    annotate(
      "text",
      x = xlim_max * 0.68,
      y = 0.98,
      label = p_label,
      hjust = 0,
      size = 4
    ) +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.text = element_text(size = 12),
      legend.position = "right",
      panel.grid = element_blank()
    )
  
  p$table <- p$table +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text = element_text(size = 11),
      panel.grid = element_blank()
    )
  
  pdf(
    file.path(output_dir, paste0(output_name, "_cbio_style.pdf")),
    width = 10,
    height = 6.5
  )
  print(p)
  dev.off()
  
  png(
    file.path(output_dir, paste0(output_name, "_cbio_style.png")),
    width = 2200,
    height = 1400,
    res = 200
  )
  print(p)
  dev.off()
  
  cat("\nLog-rank p-value:", p_value, "\n")
  
  return(list(
    data = km_df,
    fit = fit,
    logrank = logrank,
    p_value = p_value,
    plot = p,
    xlim_max = xlim_max,
    administrative_censoring = administrative_censoring
  ))
}


# ============================================================
# 13. Run final analyses
# ============================================================

results <- list()

results$OS_from_sample <- plot_km_cbio_style(
  data = df %>% filter(OS_VALID_FROM_SAMPLE),
  time_col = "OS_MONTHS_FROM_SAMPLE",
  event_col = "OS_EVENT",
  title = "Overall Survival: Primary vs Metastasis",
  x_label = "Overall Survival",
  y_label = "Overall Survival",
  output_name = "KM_OS_from_sample",
  xlim_max = X_AXIS_MAX,
  break_time_by = 5,
  administrative_censoring = ADMINISTRATIVE_CENSORING,
  show_ci = SHOW_CONFIDENCE_INTERVALS
)

results$DSS_from_sample <- plot_km_cbio_style(
  data = df %>% filter(DSS_VALID_FROM_SAMPLE),
  time_col = "DSS_MONTHS_FROM_SAMPLE",
  event_col = "DSS_EVENT",
  title = "Disease-Specific Survival: Primary vs Metastasis",
  x_label = "Disease-Specific Survival",
  y_label = "Disease-Specific Survival",
  output_name = "KM_DSS_from_sample",
  xlim_max = X_AXIS_MAX,
  break_time_by = 5,
  administrative_censoring = ADMINISTRATIVE_CENSORING,
  show_ci = SHOW_CONFIDENCE_INTERVALS
)

results$PFS_from_sample <- plot_km_cbio_style(
  data = df %>% filter(PFS_VALID_FROM_SAMPLE),
  time_col = "PFS_MONTHS_FROM_SAMPLE",
  event_col = "PFS_EVENT",
  title = "Progression-Free Survival: Primary vs Metastasis",
  x_label = "Progression-Free Survival",
  y_label = "Progression-Free Survival",
  output_name = "KM_PFS_from_sample",
  xlim_max = X_AXIS_MAX,
  break_time_by = 5,
  administrative_censoring = ADMINISTRATIVE_CENSORING,
  show_ci = SHOW_CONFIDENCE_INTERVALS
)


# ============================================================
# 14. Save p-value summary
# ============================================================

extract_result <- function(name, result) {
  if (is.null(result)) {
    return(data.frame(
      analysis = name,
      p_value = NA,
      x_axis_max_months = NA,
      administrative_censoring = NA,
      n_total = NA,
      n_metastasis = NA,
      n_primary = NA
    ))
  }
  
  result_data <- result$data
  
  return(data.frame(
    analysis = name,
    p_value = result$p_value,
    x_axis_max_months = result$xlim_max,
    administrative_censoring = result$administrative_censoring,
    n_total = nrow(result_data),
    n_metastasis = sum(result_data$KM_GROUP == "(A) Metastasis"),
    n_primary = sum(result_data$KM_GROUP == "(B) Primary")
  ))
}

pvalues <- bind_rows(
  lapply(names(results), function(nm) extract_result(nm, results[[nm]]))
)

write.csv(
  pvalues,
  file = file.path(output_dir, "logrank_pvalues_summary.csv"),
  row.names = FALSE
)

