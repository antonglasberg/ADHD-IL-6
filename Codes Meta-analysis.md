library(metafor)

##############################
# Figure 2: Main analysis with log-transformed data
##############################

# Load data
originalstudier <- read.csv("studierlog.csv", sep = ",", header = TRUE)

# Calculate effect sizes (SMD/Hedges' g)
datoriginal <- escalc(measure = "SMD",
                      m1i = originalstudier$meanIL6_adhd,
                      sd1i = originalstudier$sdIL6_adhd,
                      n1i = originalstudier$n_adhd,
                      m2i = originalstudier$meanIL6_control,
                      sd2i = originalstudier$sdIL6_control,
                      n2i = originalstudier$n_control)

# Fit models
result_re <- rma.mv(yi = datoriginal$yi, 
                    V = datoriginal$vi, 
                    random = ~ 1 | StudyNumber/Group, 
                    data = originalstudier,
                    method = "REML")

# Inverse variance (fixed-effect) model
result_iv <- rma(yi = datoriginal$yi, 
                 vi = datoriginal$vi, 
                 method = "FE")

# Calculate weights
w_re <- weights(result_re, type = "rowsum")
row_sums_weights_re <- (w_re / sum(w_re)) * 100

# Create forest plot (Figure 2)
forest(result_re, 
       slab = originalstudier$Studies, 
       main = "Log-Transformed IL-6 Levels in ADHD: All studies combined",
       ilab = cbind(originalstudier$n_adhd, originalstudier$n_control, 
                    formatC(row_sums_weights_re, format="f", digits=1)),
       ilab.xpos = c(-10, -7, -4),
       cex = 0.8,
       header = FALSE)

# Add column headers
mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
mtext("RE Weight(%)", side = 3, line = -1, at = -4, cex = 0.8)
mtext("SMD [95% CI]", side = 3, line = -1, at = 15, cex = 0.8)

# Print summary results
cat("\n--- META-ANALYSIS RESULTS ---\n\n")
cat("RANDOM EFFECTS MODEL:\n")
cat(sprintf("SMD estimate: %.3f\n", result_re$b))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", result_re$ci.lb, result_re$ci.ub))
cat(sprintf("p-value: %.4f\n", result_re$pval))
cat(sprintf("Tau² (total heterogeneity): %.3f\n", result_re$tau2))
cat(sprintf("I² (heterogeneity): %.1f%%\n\n", result_re$I2))

cat("FIXED-EFFECT (INVERSE VARIANCE) MODEL:\n")
cat(sprintf("SMD estimate: %.3f\n", result_iv$b))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", result_iv$ci.lb, result_iv$ci.ub))
cat(sprintf("p-value: %.4f\n", result_iv$pval))
cat("Tau²: Not applicable in fixed-effect model\n")
cat("I²: Not applicable in fixed-effect model\n")

##############################
# Function for Figures 3-5 (subgroup analyses with log-transformed data)
##############################

perform_meta_analysis <- function(file_path, plot_title) {
  # Load data
  orginalstudier <- read.csv(file_path, sep = ",", header = TRUE)
  
  # Log-transform data
  orginalstudier$logtransformmean_adhd <- log(orginalstudier$meanIL6_adhd) - 
    0.5 * log((orginalstudier$sdIL6_adhd^2 / orginalstudier$meanIL6_adhd^2) + 1)
  orginalstudier$logtransformSD_adhd <- sqrt(log(orginalstudier$sdIL6_adhd^2 / orginalstudier$meanIL6_adhd^2 + 1))
  orginalstudier$logtransformmean_control <- log(orginalstudier$meanIL6_control) - 
    0.5 * log((orginalstudier$sdIL6_control^2 / orginalstudier$meanIL6_control^2) + 1)
  orginalstudier$logtransformSD_control <- sqrt(log(orginalstudier$sdIL6_control^2 / orginalstudier$meanIL6_control^2 + 1))
  
  # Calculate effect sizes and variance
  logdat <- escalc(measure = "SMD",
                   m1i = orginalstudier$logtransformmean_adhd,
                   sd1i = orginalstudier$logtransformSD_adhd,
                   n1i = orginalstudier$n_adhd,
                   m2i = orginalstudier$logtransformmean_control,
                   sd2i = orginalstudier$logtransformSD_control,
                   n2i = orginalstudier$n_control)
  
  # Run both models
  result_re <- rma(yi = logdat$yi, vi = logdat$vi, method = "REML", data = orginalstudier)
  result_fe <- rma(yi = logdat$yi, vi = logdat$vi, method = "FE", data = orginalstudier)
  
  # Calculate weights for both models
  orginalstudier$weight_re <- (weights(result_re) / sum(weights(result_re))) * 100
  orginalstudier$weight_fe <- (weights(result_fe) / sum(weights(result_fe))) * 100
  orginalstudier$weight_re_formatted <- formatC(orginalstudier$weight_re, format = "f", digits = 1)
  orginalstudier$weight_fe_formatted <- formatC(orginalstudier$weight_fe, format = "f", digits = 1)
  
  # Create forest plot
  forest(result_re,
         slab = orginalstudier$Studies,
         ilab = cbind(orginalstudier$n_adhd,
                      orginalstudier$n_control,
                      orginalstudier$weight_re_formatted,
                      orginalstudier$weight_fe_formatted),
         ilab.xpos = c(-10, -7, -4, -1),
         main = plot_title,
         cex = 0.8)
  
  # Add column headers
  mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
  mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
  mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
  mtext("RE Weight (%)", side = 3, line = -1, at = -4, cex = 0.8)
  mtext("FE Weight (%)", side = 3, line = -1, at = -1, cex = 0.8)
  mtext("SMD [CI]", side = 3, line = -1, at = 15, cex = 0.8)
  
  # Add both summary estimates
  addpoly(result_re, row = -1, cex = 0.8, mlab = "Random-Effects Model", col = "blue")
  addpoly(result_fe, row = -2, cex = 0.8, mlab = "Fixed-Effect Model", col = "red")
  
  # Show results for both models
  cat("\n--- RANDOM-EFFECTS MODEL (REML) ---\n")
  print(result_re)
  cat("\n--- FIXED-EFFECT MODEL (FE) ---\n")
  print(result_fe)
  
  return(list(random_effects = result_re, fixed_effect = result_fe))
}

# Figure 3: Children subgroup
results_children <- perform_meta_analysis("children.csv", "Log-Transformed IL-6 Levels in ADHD: Participants Under 18 Years")

# Figure 4: Unmedicated subgroup
results_unmedicated <- perform_meta_analysis("unmedicated.csv", "Log-Transformed IL-6 Levels in ADHD: Unmedicated Individuals")

# Figure 5: Morning measurements subgroup
results_morning <- perform_meta_analysis("morning.csv", "Log-Transformed IL-6 Levels in ADHD: Morning Measurements")

##############################
# Figure 6: Main analysis with raw data
##############################

# Load data
originalstudier <- read.csv("studier.csv", sep = ",", header = TRUE)

# Calculate effect sizes (SMD/Hedges' g)
datoriginal <- escalc(measure = "SMD",
                      m1i = originalstudier$meanIL6_adhd,
                      sd1i = originalstudier$sdIL6_adhd,
                      n1i = originalstudier$n_adhd,
                      m2i = originalstudier$meanIL6_control,
                      sd2i = originalstudier$sdIL6_control,
                      n2i = originalstudier$n_control)

# Fit models
result_re <- rma.mv(yi = datoriginal$yi, 
                    V = datoriginal$vi, 
                    random = ~ 1 | StudyNumber/Group, 
                    data = originalstudier,
                    method = "REML")

# Inverse variance (fixed-effect) model
result_iv <- rma(yi = datoriginal$yi, 
                 vi = datoriginal$vi, 
                 method = "FE")

# Calculate weights
w_re <- weights(result_re, type = "rowsum")
row_sums_weights_re <- (w_re / sum(w_re)) * 100
w_iv <- weights(result_iv)
row_sums_weights_iv <- (w_iv / sum(w_iv)) * 100

# Create forest plot (Figure 6)
forest(result_re, 
       slab = originalstudier$Studies, 
       main = "IL-6 Levels in ADHD: Combined Analysis of All Studies",
       ilab = cbind(originalstudier$n_adhd, originalstudier$n_control, 
                    formatC(row_sums_weights_re, format="f", digits=1),
                    formatC(row_sums_weights_iv, format="f", digits=1)),
       ilab.xpos = c(-10, -7, -4, -1),
       cex = 0.8,
       header = FALSE)

# Add column headers
mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
mtext("RE Weight(%)", side = 3, line = -1, at = -4, cex = 0.8)
mtext("IV Weight(%)", side = 3, line = -1, at = -1, cex = 0.8)
mtext("SMD [95% CI]", side = 3, line = -1, at = 15, cex = 0.8)

# Print summary results
cat("\n--- META-ANALYSIS RESULTS ---\n\n")
cat("RANDOM EFFECTS MODEL:\n")
cat(sprintf("SMD estimate: %.3f\n", result_re$b))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", result_re$ci.lb, result_re$ci.ub))
cat(sprintf("p-value: %.4f\n", result_re$pval))
cat(sprintf("Tau² (total heterogeneity): %.3f\n", result_re$tau2))
cat(sprintf("I² (heterogeneity): %.1f%%\n\n", result_re$I2))

cat("FIXED-EFFECT (INVERSE VARIANCE) MODEL:\n")
cat(sprintf("SMD estimate: %.3f\n", result_iv$b))
cat(sprintf("95%% CI: [%.3f, %.3f]\n", result_iv$ci.lb, result_iv$ci.ub))
cat(sprintf("p-value: %.4f\n", result_iv$pval))
cat("Tau²: Not applicable in fixed-effect model\n")
cat("I²: Not applicable in fixed-effect model\n")

##############################
# Figure 7: Children subgroup with raw data
##############################

# Load data
originalstudier <- read.csv("children.csv", sep = ",", header = TRUE)

# Calculate effect sizes and variance
datoriginal <- escalc(measure = "SMD", 
                      m1i = originalstudier$meanIL6_adhd, 
                      sd1i = originalstudier$sdIL6_adhd, 
                      n1i = originalstudier$n_adhd, 
                      m2i = originalstudier$meanIL6_control, 
                      sd2i = originalstudier$sdIL6_control, 
                      n2i = originalstudier$n_control)

# Run both Random-Effects (REML) and Fixed-Effect (FE) models
result_re <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "REML", data = originalstudier)
result_fe <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "FE", data = originalstudier)

# Calculate weights for both models
originalstudier$weight_re <- (weights(result_re) / sum(weights(result_re))) * 100
originalstudier$weight_fe <- (weights(result_fe) / sum(weights(result_fe))) * 100
originalstudier$weight_re_formatted <- formatC(originalstudier$weight_re, format = "f", digits = 1)
originalstudier$weight_fe_formatted <- formatC(originalstudier$weight_fe, format = "f", digits = 1)

# Create forest plot (Figure 7)
forest(result_re,
       slab = originalstudier$Studies,
       ilab = cbind(originalstudier$n_adhd, 
                    originalstudier$n_control, 
                    originalstudier$weight_re_formatted,
                    originalstudier$weight_fe_formatted),
       ilab.xpos = c(-10, -7, -4, -1),
       main = "IL-6 Levels in ADHD: Participants Under 18 Years",
       cex = 0.8)

# Add column headers
mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
mtext("RE Weight (%)", side = 3, line = -1, at = -4, cex = 0.8)
mtext("FE Weight (%)", side = 3, line = -1, at = -1, cex = 0.8)
mtext("SMD [CI]", side = 3, line = -1, at = 14, cex = 0.8)

# Add both summary estimates
addpoly(result_fe, row = -1, cex = 0.8, mlab = "Fixed-Effect Model", annotate = FALSE, col = "blue")
addpoly(result_re, row = -2, cex = 0.8, mlab = "Random-Effects Model", annotate = FALSE, col = "red")

# Show results
cat("--- RANDOM-EFFECTS MODEL (REML) ---\n")
print(result_re)
cat("\n--- FIXED-EFFECT MODEL (IV) ---\n")
print(result_fe)

##############################
# Figure 8: Unmedicated subgroup with raw data
##############################

# Load data
originalstudier <- read.csv("unmedicated.csv", sep = ",", header = TRUE)

# Calculate effect sizes and variance
datoriginal <- escalc(measure = "SMD", 
                      m1i = originalstudier$meanIL6_adhd, 
                      sd1i = originalstudier$sdIL6_adhd, 
                      n1i = originalstudier$n_adhd, 
                      m2i = originalstudier$meanIL6_control, 
                      sd2i = originalstudier$sdIL6_control, 
                      n2i = originalstudier$n_control)

# Run both models
result_re <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "REML", data = originalstudier)
result_fe <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "FE", data = originalstudier)

# Calculate and format weights for both models
originalstudier$weight_re <- (weights(result_re) / sum(weights(result_re))) * 100
originalstudier$weight_fe <- (weights(result_fe) / sum(weights(result_fe))) * 100
originalstudier$weight_re_formatted <- formatC(originalstudier$weight_re, format = "f", digits = 1)
originalstudier$weight_fe_formatted <- formatC(originalstudier$weight_fe, format = "f", digits = 1)

# Create forest plot (Figure 8)
forest(result_re,
       slab = originalstudier$Studies,
       ilab = cbind(originalstudier$n_adhd, 
                    originalstudier$n_control,
                    originalstudier$weight_re_formatted,
                    originalstudier$weight_fe_formatted),
       ilab.xpos = c(-10, -7, -4, -1),
       main = "IL-6 Levels in ADHD: Unmedicated Individuals",
       cex = 0.8)

# Add column labels
mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
mtext("RE Weight (%)", side = 3, line = -1, at = -4, cex = 0.8)
mtext("FE Weight (%)", side = 3, line = -1, at = -1, cex = 0.8)
mtext("SMD [CI]", side = 3, line = -1, at = 14, cex = 0.8)

# Add both summary estimates
addpoly(result_fe, row = -1, cex = 0.8, mlab = "Fixed-Effect Model", annotate = FALSE, col = "blue")
addpoly(result_re, row = -2, cex = 0.8, mlab = "Random-Effects Model", annotate = FALSE, col = "red")

# Show results
cat("--- RANDOM-EFFECTS MODEL (REML) ---\n")
print(result_re)
cat("\n--- FIXED-EFFECT MODEL (IV) ---\n")
print(result_fe)

##############################
# Figure 9: Morning measurements subgroup with raw data
##############################

# Load data
originalstudier <- read.csv("morning.csv", sep = ",", header = TRUE)

# Calculate effect sizes and variance
datoriginal <- escalc(measure = "SMD", 
                      m1i = originalstudier$meanIL6_adhd, 
                      sd1i = originalstudier$sdIL6_adhd, 
                      n1i = originalstudier$n_adhd, 
                      m2i = originalstudier$meanIL6_control, 
                      sd2i = originalstudier$sdIL6_control, 
                      n2i = originalstudier$n_control)

# Run both models
result_re <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "REML", data = originalstudier)
result_fe <- rma(yi = datoriginal$yi, vi = datoriginal$vi, method = "FE", data = originalstudier)

# Calculate and format weights
originalstudier$weight_re <- (weights(result_re) / sum(weights(result_re))) * 100
originalstudier$weight_fe <- (weights(result_fe) / sum(weights(result_fe))) * 100
originalstudier$weight_re_formatted <- formatC(originalstudier$weight_re, format = "f", digits = 1)
originalstudier$weight_fe_formatted <- formatC(originalstudier$weight_fe, format = "f", digits = 1)

# Create forest plot (Figure 9)
forest(result_re,
       slab = originalstudier$Studies,
       ilab = cbind(originalstudier$n_adhd, 
                    originalstudier$n_control,
                    originalstudier$weight_re_formatted,
                    originalstudier$weight_fe_formatted),
       ilab.xpos = c(-10, -7, -4, -1),
       main = "IL-6 Levels in ADHD: Morning Measurements",
       cex = 0.8)

# Add column labels
mtext("Author and Year", side = 3, line = -1, at = -18, cex = 0.8)
mtext("ADHD(n)", side = 3, line = -1, at = -10, cex = 0.8)
mtext("Control(n)", side = 3, line = -1, at = -7, cex = 0.8)
mtext("RE Weight (%)", side = 3, line = -1, at = -4, cex = 0.8)
mtext("FE Weight (%)", side = 3, line = -1, at = -1, cex = 0.8)
mtext("SMD [CI]", side = 3, line = -1, at = 14, cex = 0.8)

# Add both summary estimates
addpoly(result_fe, row = -1, cex = 0.8, mlab = "Fixed-Effect Model", annotate = FALSE, col = "blue")
addpoly(result_re, row = -2, cex = 0.8, mlab = "Random-Effects Model", annotate = FALSE, col = "red")

# Show results
cat("--- RANDOM-EFFECTS MODEL (REML) ---\n")
print(result_re)
cat("\n--- FIXED-EFFECT MODEL (IV) ---\n")
print(result_fe)
