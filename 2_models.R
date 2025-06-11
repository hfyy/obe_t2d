install.packages("msm")
install.packages("elect")

library(msm)
library(elect)
library(haven)
library(dplyr)
library(ggplot2)
library(tibble)

data <- read_dta("/Users/yiyue/Library/CloudStorage/OneDrive-CUHK-Shenzhen/T2D_obesity/data/sample_t2dobe.dta")
data <- data %>% select(c("wave","raeduc","hhidpn","rbmi","rdiab_raw","race","age","gender","rdiab","robe","dead","hhidpn",
                          "PGS_BMI","age_base","ec2_std","PGS_T2D","wave","year"))
data <- data %>% mutate(age = age-45)
data <- data %>% mutate(edu = case_when(raeduc == 1 ~ 0,
                                        raeduc == 2 ~ 0, 
                                        raeduc == 3 ~1 ))

data <- data %>% mutate(state = case_when(rdiab_raw ==0 & robe ==0 & dead==0~ 1, ## (1) no obese or T2D
                                          rdiab_raw ==1 & robe==0  & dead==0  ~2,    ## (2) T2D, no obese
                                          rdiab_raw ==0 & robe==1  & dead==0 ~ 3,   ## (3) obese, no T2D
                                          rdiab_raw ==1 & robe==1  & dead==0 ~ 4, ## (4) T2D, obese    
                                          dead ==1 ~ 5))  %>%   
                 arrange(hhidpn, wave)


statetable.msm(state, hhidpn, data=data)

any(!is.finite(data$state))
any(!is.finite(data$year))

data <- data[is.finite(data$state) , ]


#################################################################################
##               baseline
#################################################################################

Q  <-  rbind ( c(0, 0.02,0.06,0,0.02),
               c(0.04, 0, 0,0.05,0.04),
               c(0.2,0, 0, 0.06,0.02),
               c(0,0.12,0.05, 0,0.05),
               c(0,0,0,0,0))

msm <- msm( state ~ age, subject=hhidpn, data = data, death = 5, qmatrix = Q, 
            control=list(fnscale=10000,maxit=5000),  center=FALSE)
   msm
 

#################################################################################
##                model with covariates and interactions                               
################################################################################# 
 
msm2 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
               covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D+ ec2_std + factor(edu) )
   msm2

     
  pmatrix <- data.frame(pmatrix.msm(msm2, covariates = list(age=0))) 
  
   pmatrix.msm(msm2, covariates = list(age=2))
   pmatrix.msm(msm2, covariates = list(age=4))
   pmatrix.msm(msm2, covariates = list(age=6))
   pmatrix.msm(msm2, covariates = list(age=8))
   pmatrix.msm(msm2, covariates = list(age=10))
   pmatrix.msm(msm2, covariates = list(age=12))
   pmatrix.msm(msm2, covariates = list(age=14))
   pmatrix.msm(msm2, covariates = list(age=16))
   pmatrix.msm(msm2, covariates = list(age=18))
   pmatrix.msm(msm2, covariates = list(age=20))
   pmatrix.msm(msm2, covariates = list(age=22))
   pmatrix.msm(msm2, covariates = list(age=24))
   pmatrix.msm(msm2, covariates = list(age=26))
   pmatrix.msm(msm2, covariates = list(age=28))
   pmatrix.msm(msm2, covariates = list(age=30))
   pmatrix.msm(msm2, covariates = list(age=32))
   pmatrix.msm(msm2, covariates = list(age=34))
   pmatrix.msm(msm2, covariates = list(age=36))
   pmatrix.msm(msm2, covariates = list(age=38))
   pmatrix.msm(msm2, covariates = list(age=40))
              
  
msm2_inter <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
                   covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D + ec2_std + ec2_std*PGS_BMI + factor(edu) )

msm2_inter2 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, death=5, control=list(fnscale=5000,maxit=500), 
                   covariates= ~ age + factor(gender) + PGS_BMI +PGS_T2D + ec2_std + ec2_std*PGS_T2D + factor(edu) )



#################################################################################
##               hidden Markov model                  
################################################################################# 

Q  <-  rbind ( c(0, 0.01,0.03,0,0.1),
               c(0, 0, 0,0.05,0.04),
               c(0.1,0, 0, 0.03,0.01),
               c(0,0.06,0, 0,0.02),
               c(0,0,0,0,0))

ematrix <- rbind(c(0, 0.001, 0.001, 0.001,0),
                 c(0.5, 0, 0.001, 0.001,0),
                 c(0.001, 0.5, 0, 0.001,0),
                 c(0.5, 0.001, 0.5, 0,0),
                 c(0,0,0,0,0))


msm3 <- msm( state ~ time, subject=hhidpn, data = data, death = 5, qmatrix = Q, ematrix = ematrix,covariates = ~age, 
            control=list(fnscale=10000,maxit=5000,reltol= 1e-08 ),  center=FALSE, method = "BFGS")
 msm3

msm4 <- msm( state ~ time, subject=hhidpn, data = data, qmatrix = Q, ematrix = ematrix, death=5, control=list(fnscale=5000,maxit=500, reltol = 1e-08),method = "BFGS", 
             covariates= ~ age + factor(gender) + PGS_BMI + PGS_T2D + ec2_std + factor(edu))
 msm4
  
 #################################################################################
 ##             to plot the figures                 
 ################################################################################# 
 #### Figure 2, baseline 
 # Load necessary library
 library(ggplot2)
 
 # Create the data frame
 data <- data.frame(
   Transition = c(
     "nO/nD-nO/nD", "nO/nD-nO/D", "nO/nD-O/nD", "nO/nD-death",
     "nO/D-nO/nD", "nO/D-nO/D", "nO/D-O/D", "nO/D-death",
     "O/nD-nO/nD", "O/nD-O/nD", "O/nD-O/D", "O/nD-death",
     "O/D-nO/D", "O/D-O/nD", "O/D-O/D", "O/D-death"
   ),
   Baseline = c(
     0.052, 0.010, 0.029, 0.012,
     0.021, 0.107, 0.048, 0.038,
     0.075, 0.115, 0.031, 0.008,
     0.068, 0.013, 0.102, 0.021
   ),
   SE = c(
     0.001, 0.000, 0.001, 0.000,
     0.001, 0.003, 0.002, 0.002,
     0.002, 0.002, 0.001, 0.001,
     0.003, 0.001, 0.003, 0.001
   )
 )
 
 # Calculate the Confidence Intervals
 data$Lower_CI <- data$Baseline - 1.96 * data$SE
 data$Upper_CI <- data$Baseline + 1.96 * data$SE
 
 # Add the Group variable
 data$Group <- c(
   rep("nO/nD", 4),
   rep("nO/D", 4),
   rep("O/nD", 4),
   rep("O/D", 4)
 )
 
 # Ensure the Transition variable is ordered as in the original vector
 data$Transition <- factor(data$Transition, levels = rev(c(
   "nO/nD-nO/nD", "nO/nD-nO/D", "nO/nD-O/nD", "nO/nD-death",
   "nO/D-nO/nD", "nO/D-nO/D", "nO/D-O/D", "nO/D-death",
   "O/nD-nO/nD", "O/nD-O/nD", "O/nD-O/D", "O/nD-death",
   "O/D-nO/D", "O/D-O/nD", "O/D-O/D", "O/D-death"
 )))
 
 # Ensure the Group variable is ordered as desired
 data$Group <- factor(data$Group, levels = c("nO/nD", "nO/D", "O/nD", "O/D"))
 
 # Create the horizontal bar plot
 ggplot(data, aes(x = Baseline, y = Transition, fill = Baseline > 0)) +
   # Horizontal bars
   geom_bar(stat = "identity", color = "black") +
   # Error bars
   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, color = "red") +
   # Separate groups with facet_grid
   facet_grid(Group ~ ., scales = "free", space = "free") +
   # Customize the color palette for bars
   scale_fill_manual(values = c("skyblue", "skyblue"), guide = "none") +
   # Add labels and customize themes
   labs(
     title = "",
     x = "Baseline hazard rates",
     y = "Transitions"
   ) +
   theme_minimal() +
   theme(
     strip.text.y.right = element_blank(),  # Remove facet strip text
     panel.grid.major.y = element_blank(),  # Remove default gridlines
     panel.spacing = unit(1, "lines"),      # Add spacing between groups
     axis.text.y = element_text(size = 10), # Set y-axis text size
     text = element_text(size = 14)         # Set overall text size
   )
 
 #### Figure 3, baseline hazard rates and coefficients 
 library(dplyr)
 library(tidyr)
 library(forcats)
 library(gridExtra)
 
 # Create data frame with transition matrix data
 transitions <- c(
   "nO/nD-nO/nD", "nO/nD-nO/D", "nO/nD-O/nD", "nO/nD-death",
   "nO/D-nO/nD", "nO/D-nO/D", "nO/D-O/D", "nO/D-death",
   "O/nD-nO/nD", "O/nD-O/nD", "O/nD-O/D", "O/nD-death",
   "O/D-nO/D", "O/D-O/nD", "O/D-O/D", "O/D-death"
 )
 
 # Create data frame with coefficients and standard errors
 transition_data <- data.frame(
   transition = rep(transitions, each = 7),
   variable = rep(c("Baseline", "Age", "Female", "PRS_BMI", "PRS_T2D", "Early_Condition", "Bachelors_Degree"), length(transitions)),
   coefficient = c(
     # nO/nD-nO/nD
     -0.042, NA, NA, NA, NA, NA, NA,
     # nO/nD-nO/D
     0.009, 0.009, -0.201, -0.004, 0.071, -0.013, -0.088,
     # nO/nD-O/nD
     0.029, -0.011, -0.03, 0.154, -0.033, -0.013, -0.051,
     # nO/nD-death
     0.005, 0.058, -0.188, 0.056, -0.051, -0.022, -0.080,
     # nO/D-nO/nD
     0.023, -0.011, 0.105, -0.066, -0.02, -0.051, 0.201,
     # nO/D-nO/D
     -0.087, NA, NA, NA, NA, NA, NA,
     # nO/D-O/D
     0.050, -0.011, -0.007, 0.100, 0.038, 0.081, -0.150,
     # nO/D-death
     0.014, 0.046, -0.09, 0.064, -0.029, -0.008, -0.059,
     # O/nD-nO/nD
     0.081, 0.009, -0.011, -0.056, 0.014, 0.038, 0.031,
     # O/nD-O/nD
     -0.116, NA, NA, NA, NA, NA, NA,
     # O/nD-O/D
     0.030, -0.001, -0.102, 0, 0.055, -0.025, -0.043,
     # O/nD-death
     0.005, 0.057, -0.087, 0.006, -0.023, -0.046, -0.025,
     # O/D-nO/D
     0.069, 0.013, -0.025, -0.097, 0.025, 0.033, -0.081,
     # O/D-O/nD
     0.014, -0.014, 0.054, -0.025, -0.078, -0.031, 0.032,
     # O/D-O/D
     -0.097, NA, NA, NA, NA, NA, NA,
     # O/D-death
     0.014, 0.039, -0.170, 0.007, -0.005, 0.035, -0.122
   ),
   se = c(
     # nO/nD-nO/nD
     0.001, NA, NA, NA, NA, NA, NA,
     # nO/nD-nO/D
     0.000, 0.002, 0.031, 0.016, 0.015, 0.019, 0.034,
     # nO/nD-O/nD
     0.021, 0.001, 0.018, 0.009, 0.009, 0.011, 0.019,
     # nO/nD-death
     0.000, 0.002, 0.028, 0.014, 0.014, 0.018, 0.030,
     # nO/D-nO/nD
     0.002, 0.003, 0.060, 0.031, 0.033, 0.037, 0.065,
     # nO/D-nO/D
     0.003, NA, NA, NA, NA, NA, NA,
     # nO/D-O/D
     0.003, 0.002, 0.042, 0.023, 0.021, 0.026, 0.047,
     # nO/D-death
     0.002, 0.003, 0.051, 0.026, 0.026, 0.035, 0.053,
     # O/nD-nO/nD
     0.002, 0.001, 0.020, 0.011, 0.011, 0.013, 0.022,
     # O/nD-O/nD
     0.002, NA, NA, NA, NA, NA, NA,
     # O/nD-O/D
     0.001, 0.002, 0.031, 0.017, 0.015, 0.019, 0.033,
     # O/nD-death
     0.001, 0.004, 0.070, 0.038, 0.031, 0.043, 0.076,
     # O/D-nO/D
     0.003, 0.002, 0.036, 0.020, 0.019, 0.023, 0.039,
     # O/D-O/nD
     0.001, 0.005, 0.080, 0.043, 0.043, 0.050, 0.085,
     # O/D-O/D
     0.004, NA, NA, NA, NA, NA, NA,
     # O/D-death
     0.002, 0.004, 0.072, 0.038, 0.035, 0.046, 0.078
   )
 )
 
 # Calculate significance for plotting
 transition_data <- transition_data %>%
   mutate(
     signif = case_when(
       is.na(se) ~ "none",
       abs(coefficient)/se >= 3.291 ~ "***",
       abs(coefficient)/se >= 2.576 ~ "**",
       abs(coefficient)/se >= 1.96 ~ "*",
       TRUE ~ ""
     )
   )
 
 # Create a better ordering for transitions
 ordered_transitions <- c(
   "nO/nD-nO/D", "nO/nD-O/nD", "nO/nD-death",
   "nO/D-nO/nD", "nO/D-O/D", "nO/D-death",
   "O/nD-nO/nD", "O/nD-O/D", "O/nD-death",
   "O/D-nO/D", "O/D-O/nD", "O/D-death"
 )
 # Filter to include only off-diagonal transitions
 filtered_data <- transition_data %>%
   filter(transition %in% ordered_transitions)
 
 # Factor the transitions for consistent ordering
 filtered_data$transition <- factor(filtered_data$transition, 
                                    levels = rev(ordered_transitions))  # Reverse for proper top-to-bottom display
 
 # Create more readable variable names for the plot
 filtered_data$variable <- factor(filtered_data$variable, 
                                  levels = c("Baseline", "Age", "Female", "PRS_BMI", "PRS_T2D", 
                                             "Early_Condition", "Bachelors_Degree"),
                                  labels = c("Baseline", "Age", "Female", "PRS for BMI", "PRS for T2D", 
                                             "Early Condition", "Bachelor's"))
 
 # Create a multi-panel plot with different geoms for Baseline vs other variables
 p <- ggplot(filtered_data, aes(x = coefficient, y = transition)) +
   # Add reference line at zero
   geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
   
   # Use bars for Baseline
   geom_bar(data = filtered_data %>% filter(variable == "Baseline"),
            stat = "identity", 
            aes(fill = coefficient > 0)) +
   
   # Use error bars and points for other variables
   geom_errorbar(data = filtered_data %>% filter(variable != "Baseline"),
                 aes(xmin = coefficient - 1.96*se, 
                     xmax = coefficient + 1.96*se), 
                 height = 0.3, width = 0) +
   geom_point(data = filtered_data %>% filter(variable != "Baseline"),
              aes(color = coefficient > 0), 
              size = 2) +
   
   # Add significance indicators
   geom_text(aes(label = signif, 
                 x = ifelse(coefficient > 0, 
                            ifelse(variable == "Baseline", 
                                   coefficient + 0.005, 
                                   coefficient + 1.96*se + 0.01), 
                            ifelse(variable == "Baseline", 
                                   coefficient - 0.005, 
                                   coefficient - 1.96*se - 0.01))),
             hjust = ifelse(filtered_data$coefficient > 0, 0, 1),
             size = 2.5) +
   
   # Set colors for points and bars
   scale_color_manual(values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D"), guide = "none") +
   scale_fill_manual(values = c("TRUE" = "#4393C3", "FALSE" = "#D6604D"), guide = "none") +
   
   # Create facets for each variable
   facet_wrap(~ variable, scales = "free_x", ncol = 4) +
   
   # Add labels
   labs(x = "Coefficient", y = "Transition") +
   
   # Set theme
   theme_minimal() +
   theme(
     strip.background = element_rect(fill = "lightgray", color = NA),
     strip.text = element_text(face = "bold", size = 9),
     axis.text.y = element_text(size = 7),
     axis.text.x = element_text(size = 7),
     axis.title = element_text(size = 9),
     panel.grid.major.y = element_blank(),
     panel.spacing = unit(0.5, "lines")
   )
 
 # Display the plot
 print(p)
 
 #### Figure 4. 
 # Create the data frame with cleaned column names
 
 data <- data.frame(
   Covariate = c(
     "PRS for BMI and T2D", "PRS for BMI and T2D",
     "Education", "Education",
     "Gender", "Gender"
   ),
   Level = c(
     "-1 std. dev.", "+1 std. dev.",
     "< Bachelor's", ">= Bachelor's",
     "Male", "Female"
   ),
   All = c(
     25.6, 25.4,
     25.7, 25.3,
     25.2, 25.8
   ),
   nO_nD = c(
     17.1, 15.8,
     16.2, 16.8,
     16.1, 16.8
   ),
   nO_D = c(
     4.2, 4.2,
     4.5, 3.7,
     4.2, 4.1
   ),
   O_nD  = c(
     2.7, 3.2,
     3.0, 2.9,
     2.9, 3.0
   ),
   O_D = c(
     1.6, 2.2,
     2.0, 1.9,
     2.0, 1.9
   )
 )
 
 # Set the Covariate variable as a factor with a specific order
 data$Covariate <- factor(data$Covariate, levels = c("Education", "Gender", "PRS for BMI and T2D"))
 
 # Reshape the data into long format
 data_long <- data %>%
   pivot_longer(
     cols = c("nO_nD", "nO_D", "O_nD", "O_D"),  # Exclude "All"
     names_to = "State",
     values_to = "Value"
   )
 
 # Set the State variable as a factor with the desired order
 data_long$State <- factor(data_long$State, levels = c("nO_nD", "nO_D", "O_nD", "O_D"))
 
 # Create the horizontal stacked bar graph
 ggplot(data_long, aes(x = Value, y = Level, fill = State)) +
   geom_bar(stat = "identity", position = "stack", color = "black", width = 0.9) +  # Bars take up full row space
   geom_text(
     aes(label = Value), 
     position = position_stack(vjust = 0.5),  # Place labels in the middle of each segment
     size = 3,  # Adjust text size
     color = "black"  # Text color
   ) +
   facet_wrap(~ Covariate, scales = "free_y", ncol = 1) +  # One column, ordered by Covariate factor
   scale_fill_manual(
     values = c(
       "nO_nD" = "#4E79A7",  # Blue
       "nO_D" = "#F28E2B",   # Orange
       "O_nD" = "#76B7B2",   # Teal
       "O_D" = "#E15759"     # Red
     ),
     name = "State",
     labels = c("nO/nD", "nO/D", "O/nD", "O/D")  # Custom labels for clarity
   ) +
   labs(
     title = "",
     x = "Years Lived",  # Adjusted label to reflect actual years lived
     y = "Covariate Level"
   ) +
   theme_minimal() +
   theme(
     strip.text = element_text(size = 12, face = "bold"),  # Facet labels
     axis.text.y = element_text(size = 10),  # Y-axis text size
     axis.text.x = element_text(size = 10),  # X-axis text size
     legend.position = "bottom",  # Move legend to the bottom
     legend.title = element_text(size = 12),  # Legend title size
     legend.text = element_text(size = 10),  # Legend text size
     legend.key.width = unit(1, "cm"),  # Widen legend keys
     legend.key.height = unit(0.5, "cm")  # Reduce legend key height
   )