# Analyze HPLC data 
library(dplyr)
library(ggplot2)
# Load datasets
setwd('../../data/experimentalData/')

# Saccharomyces cerevisiae
filename <- 'HPLC_cerevisiae.csv'
data_cerevisiae <- read.csv(filename,header = TRUE,sep = ';',dec = ',')

# Create semi-log plots for replicates to check for exponential growth
# OD replicte 1
ggplot(data_cerevisiae,aes(x=Time,y=log(OD_repl1))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# DW replicate 1
ggplot(data_cerevisiae,aes(x=Time,y=log(CDW_repl1))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# OD replicte 2
ggplot(data_cerevisiae,aes(x=Time,y=log(OD_repl2))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# DW replicate 2
ggplot(data_cerevisiae,aes(x=Time,y=log(CDW_repl2))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# OD replicte 3
ggplot(data_cerevisiae,aes(x=Time,y=log(OD_repl3))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# DW replicate 3
ggplot(data_cerevisiae,aes(x=Time,y=log(CDW_repl3))) +
  geom_point() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  scale_x_continuous(limits = c(0,16),breaks = c(seq(1,16,by=1)))
# Based on these plots, use t3, t4, t6, t7, t8 from repl_1,
# t3, t5, t6, t7, t8 from repl_2 and t4, t5, t6, t7, t8 from repl_3

cerevisiae_repl1 <- data_cerevisiae[grep("t3|t4|t5|t6|t7|t8",data_cerevisiae$Timepoint),]
cerevisiae_repl2 <- data_cerevisiae[grep("t5|t6|t7|t8|t9",data_cerevisiae$Timepoint),]
cerevisiae_repl3 <- data_cerevisiae[grep("t4|t5|t6|t7|t8",data_cerevisiae$Timepoint),]

# Define molecular weight of metabolites
Glc_mw <- 180.16/1000 #g/mmol
EtOH_mw <- 46.07/1000 #g/mmol
Glycerol_mw <- 92.09/1000 #g/mmol
Pyr_mw <- 87.05/1000 #g/mmol
Ace_mw <- 59.04/1000 #g/mmol
Succ_mw <- 116.07/1000 #g/mmol

# growth rate based on dry weight
# Replicate 1
y_cer_1 <- log(cerevisiae_repl1$CDW_repl1)
x_cer_1 <- cerevisiae_repl1$Time
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
grate_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
grate_cer_1 <- substring(grate_cer_1,1,7)
r_squared_grate_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_grate_cer_1 <- substring(r_squared_grate_cer_1,1,6)
max_y_grate_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
max_x_grate_cer_1 <- max(cerevisiae_repl1$Time)
# Plot curve
ggplot(cerevisiae_repl1, aes(x=Time , y=log(CDW_repl1))) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_grate_cer_1*(0.7), y=max_y_grate_cer_1*(0.8), label= r_squared_grate_cer_1) + 
  annotate("text", x=max_x_grate_cer_1*(0.7), y=max_y_grate_cer_1*(0.9), label= grate_cer_1) + 
  annotate("text", x=max_x_grate_cer_1*(0.9), y=max_y_grate_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_grate_cer_1*(0.9), y=max_y_grate_cer_1*(0.9), label= 'h-1') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# Replicate 2
y_cer_2 <- log(cerevisiae_repl2$CDW_repl2)
x_cer_2 <- cerevisiae_repl2$Time
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
grate_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
grate_cer_2 <- substring(grate_cer_2,1,7)
r_squared_grate_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_grate_cer_2 <- substring(r_squared_grate_cer_2,1,6)
max_y_grate_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
max_x_grate_cer_2 <- max(cerevisiae_repl2$Time)
# Plot curve
ggplot(cerevisiae_repl2, aes(x=Time , y=log(CDW_repl2))) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_grate_cer_2*(0.7), y=max_y_grate_cer_2*(0.8), label= r_squared_grate_cer_2) + 
  annotate("text", x=max_x_grate_cer_2*(0.7), y=max_y_grate_cer_2*(0.9), label= grate_cer_2) + 
  annotate("text", x=max_x_grate_cer_2*(0.9), y=max_y_grate_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_grate_cer_2*(0.9), y=max_y_grate_cer_2*(0.9), label= 'h-1') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# Replicate 3
y_cer_3 <- log(cerevisiae_repl3$CDW_repl3)
x_cer_3 <- cerevisiae_repl3$Time
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
grate_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
grate_cer_3 <- substring(grate_cer_3,1,7)
r_squared_grate_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_grate_cer_3 <- substring(r_squared_grate_cer_3,1,6)
max_y_grate_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
max_x_grate_cer_3 <- max(cerevisiae_repl3$Time)
# Plot curve
ggplot(cerevisiae_repl3, aes(x=Time , y=log(CDW_repl3))) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_grate_cer_3*(0.7), y=max_y_grate_cer_3*(0.8), label= r_squared_grate_cer_3) + 
  annotate("text", x=max_x_grate_cer_3*(0.7), y=max_y_grate_cer_3*(0.9), label= grate_cer_3) + 
  annotate("text", x=max_x_grate_cer_3*(0.9), y=max_y_grate_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_grate_cer_3*(0.9), y=max_y_grate_cer_3*(0.9), label= 'h-1') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# glucose uptake rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Glucose_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
glucose_uptake_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
glucose_uptake_cer_1 <- substring(glucose_uptake_cer_1,1,7)
r_squared_glucose_uptake_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_glucose_uptake_cer_1 <- substring(r_squared_glucose_uptake_cer_1,1,6)
max_y_glucose_uptake_cer_1 <- max(cerevisiae_repl1$Glucose_repl1)
max_x_glucose_uptake_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
glucose_uptake_rate_cer_1 <- (as.numeric(glucose_uptake_cer_1)/Glc_mw)*as.numeric(grate_cer_1)

# plot glucose uptake yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1 , y=Glucose_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glucose_uptake_cer_1*(0.7), y=max_y_glucose_uptake_cer_1*(0.8), label= r_squared_glucose_uptake_cer_1) + 
  annotate("text", x=max_x_glucose_uptake_cer_1*(0.7), y=max_y_glucose_uptake_cer_1*(0.9), label= glucose_uptake_cer_1) + 
  annotate("text", x=max_x_glucose_uptake_cer_1*(0.9), y=max_y_glucose_uptake_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glucose_uptake_cer_1*(0.9), y=max_y_glucose_uptake_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Glucose_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
glucose_uptake_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
glucose_uptake_cer_2 <- substring(glucose_uptake_cer_2,1,7)
r_squared_glucose_uptake_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_glucose_uptake_cer_2 <- substring(r_squared_glucose_uptake_cer_2,1,6)
max_y_glucose_uptake_cer_2 <- max(cerevisiae_repl2$Glucose_repl2)
max_x_glucose_uptake_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
glucose_uptake_rate_cer_2 <- (as.numeric(glucose_uptake_cer_2)/Glc_mw)*as.numeric(grate_cer_2)

# plot glucose uptake yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2 , y=Glucose_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glucose_uptake_cer_2*(0.7), y=max_y_glucose_uptake_cer_2*(0.8), label= r_squared_glucose_uptake_cer_2) + 
  annotate("text", x=max_x_glucose_uptake_cer_2*(0.7), y=max_y_glucose_uptake_cer_2*(0.9), label= glucose_uptake_cer_2) + 
  annotate("text", x=max_x_glucose_uptake_cer_2*(0.9), y=max_y_glucose_uptake_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glucose_uptake_cer_2*(0.9), y=max_y_glucose_uptake_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Glucose_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
glucose_uptake_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
glucose_uptake_cer_3 <- substring(glucose_uptake_cer_3,1,7)
r_squared_glucose_uptake_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_glucose_uptake_cer_3 <- substring(r_squared_glucose_uptake_cer_3,1,6)
max_y_glucose_uptake_cer_3 <- max(cerevisiae_repl3$Glucose_repl3)
max_x_glucose_uptake_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
glucose_uptake_rate_cer_3 <- (as.numeric(glucose_uptake_cer_3)/Glc_mw)*as.numeric(grate_cer_3)

# plot glucose uptake yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3 , y=Glucose_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glucose_uptake_cer_3*(0.7), y=max_y_glucose_uptake_cer_3*(0.8), label= r_squared_glucose_uptake_cer_3) + 
  annotate("text", x=max_x_glucose_uptake_cer_3*(0.7), y=max_y_glucose_uptake_cer_3*(0.9), label= glucose_uptake_cer_3) + 
  annotate("text", x=max_x_glucose_uptake_cer_3*(0.9), y=max_y_glucose_uptake_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glucose_uptake_cer_3*(0.9), y=max_y_glucose_uptake_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# ethanol secretion rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Ethanol_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
ethanol_secretion_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
ethanol_secretion_cer_1 <- substring(ethanol_secretion_cer_1,1,7)
r_squared_ethanol_secretion_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_ethanol_secretion_cer_1 <- substring(r_squared_ethanol_secretion_cer_1,1,6)
max_y_ethanol_secretion_cer_1 <- max(cerevisiae_repl1$Ethanol_repl1)
max_x_ethanol_secretion_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
ethanol_secretion_rate_cer_1 <- (as.numeric(ethanol_secretion_cer_1)/EtOH_mw)*as.numeric(grate_cer_1)

# plot ethanol secretion yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1 , y=Ethanol_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_ethanol_secretion_cer_1*(0.7), y=max_y_ethanol_secretion_cer_1*(0.8), label= r_squared_ethanol_secretion_cer_1) + 
  annotate("text", x=max_x_ethanol_secretion_cer_1*(0.7), y=max_y_ethanol_secretion_cer_1*(0.9), label= ethanol_secretion_cer_1) + 
  annotate("text", x=max_x_ethanol_secretion_cer_1*(0.9), y=max_y_ethanol_secretion_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_ethanol_secretion_cer_1*(0.9), y=max_y_ethanol_secretion_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Ethanol_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
ethanol_secretion_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
ethanol_secretion_cer_2 <- substring(ethanol_secretion_cer_2,1,7)
r_squared_ethanol_secretion_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_ethanol_secretion_cer_2 <- substring(r_squared_ethanol_secretion_cer_2,1,6)
max_y_ethanol_secretion_cer_2 <- max(cerevisiae_repl2$Ethanol_repl2)
max_x_ethanol_secretion_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
ethanol_secretion_rate_cer_2 <- (as.numeric(ethanol_secretion_cer_2)/EtOH_mw)*as.numeric(grate_cer_2)

# plot ethanol secretion yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2 , y=Ethanol_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_ethanol_secretion_cer_2*(0.7), y=max_y_ethanol_secretion_cer_2*(0.8), label= r_squared_ethanol_secretion_cer_2) + 
  annotate("text", x=max_x_ethanol_secretion_cer_2*(0.7), y=max_y_ethanol_secretion_cer_2*(0.9), label= ethanol_secretion_cer_2) + 
  annotate("text", x=max_x_ethanol_secretion_cer_2*(0.9), y=max_y_ethanol_secretion_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_ethanol_secretion_cer_2*(0.9), y=max_y_ethanol_secretion_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Ethanol_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
ethanol_secretion_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
ethanol_secretion_cer_3 <- substring(ethanol_secretion_cer_3,1,7)
r_squared_ethanol_secretion_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_ethanol_secretion_cer_3 <- substring(r_squared_ethanol_secretion_cer_3,1,6)
max_y_ethanol_secretion_cer_3 <- max(cerevisiae_repl3$Ethanol_repl3)
max_x_ethanol_secretion_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
ethanol_secretion_rate_cer_3 <- (as.numeric(ethanol_secretion_cer_3)/EtOH_mw)*as.numeric(grate_cer_3)

# plot ethanol secretion yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3 , y=Ethanol_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_ethanol_secretion_cer_3*(0.7), y=max_y_ethanol_secretion_cer_3*(0.8), label= r_squared_ethanol_secretion_cer_3) + 
  annotate("text", x=max_x_ethanol_secretion_cer_3*(0.7), y=max_y_ethanol_secretion_cer_3*(0.9), label= ethanol_secretion_cer_3) + 
  annotate("text", x=max_x_ethanol_secretion_cer_3*(0.9), y=max_y_ethanol_secretion_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_ethanol_secretion_cer_3*(0.9), y=max_y_ethanol_secretion_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# glycerol secretion rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Glycerol_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
glycerol_secretion_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
glycerol_secretion_cer_1 <- substring(glycerol_secretion_cer_1,1,7)
r_squared_glycerol_secretion_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_glycerol_secretion_cer_1 <- substring(r_squared_glycerol_secretion_cer_1,1,6)
max_y_glycerol_secretion_cer_1 <- max(cerevisiae_repl1$Glycerol_repl1)
max_x_glycerol_secretion_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
glycerol_secretion_rate_cer_1 <- (as.numeric(glycerol_secretion_cer_1)/Glycerol_mw)*as.numeric(grate_cer_1)

# plot glycerol secretion yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1 , y=Glycerol_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glycerol_secretion_cer_1*(0.7), y=max_y_glycerol_secretion_cer_1*(0.8), label= r_squared_glycerol_secretion_cer_1) + 
  annotate("text", x=max_x_glycerol_secretion_cer_1*(0.7), y=max_y_glycerol_secretion_cer_1*(0.9), label= glycerol_secretion_cer_1) + 
  annotate("text", x=max_x_glycerol_secretion_cer_1*(0.9), y=max_y_glycerol_secretion_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glycerol_secretion_cer_1*(0.9), y=max_y_glycerol_secretion_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Glycerol_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
glycerol_secretion_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
glycerol_secretion_cer_2 <- substring(glycerol_secretion_cer_2,1,7)
r_squared_glycerol_secretion_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_glycerol_secretion_cer_2 <- substring(r_squared_glycerol_secretion_cer_2,1,6)
max_y_glycerol_secretion_cer_2 <- max(cerevisiae_repl2$Glycerol_repl2)
max_x_glycerol_secretion_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
glycerol_secretion_rate_cer_2 <- (as.numeric(glycerol_secretion_cer_2)/Glycerol_mw)*as.numeric(grate_cer_2)

# plot glycerol secretion yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2, y=Glycerol_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glycerol_secretion_cer_2*(0.7), y=max_y_glycerol_secretion_cer_2*(0.8), label= r_squared_glycerol_secretion_cer_2) + 
  annotate("text", x=max_x_glycerol_secretion_cer_2*(0.7), y=max_y_glycerol_secretion_cer_2*(0.9), label= glycerol_secretion_cer_2) + 
  annotate("text", x=max_x_glycerol_secretion_cer_2*(0.9), y=max_y_glycerol_secretion_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glycerol_secretion_cer_2*(0.9), y=max_y_glycerol_secretion_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Glycerol_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
glycerol_secretion_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
glycerol_secretion_cer_3 <- substring(glycerol_secretion_cer_3,1,7)
r_squared_glycerol_secretion_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_glycerol_secretion_cer_3 <- substring(r_squared_glycerol_secretion_cer_3,1,6)
max_y_glycerol_secretion_cer_3 <- max(cerevisiae_repl3$Glycerol_repl3)
max_x_glycerol_secretion_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
glycerol_secretion_rate_cer_3 <- (as.numeric(glycerol_secretion_cer_3)/Glycerol_mw)*as.numeric(grate_cer_3)

# plot glycerol secretion yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3, y=Glycerol_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_glycerol_secretion_cer_3*(0.7), y=max_y_glycerol_secretion_cer_3*(0.8), label= r_squared_glycerol_secretion_cer_3) + 
  annotate("text", x=max_x_glycerol_secretion_cer_3*(0.7), y=max_y_glycerol_secretion_cer_3*(0.9), label= glycerol_secretion_cer_3) + 
  annotate("text", x=max_x_glycerol_secretion_cer_3*(0.9), y=max_y_glycerol_secretion_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_glycerol_secretion_cer_3*(0.9), y=max_y_glycerol_secretion_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# pyruvate secretion rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Pyruvate_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
pyruvate_secretion_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
pyruvate_secretion_cer_1 <- substring(pyruvate_secretion_cer_1,1,7)
r_squared_pyruvate_secretion_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_pyruvate_secretion_cer_1 <- substring(r_squared_pyruvate_secretion_cer_1,1,6)
max_y_pyruvate_secretion_cer_1 <- max(cerevisiae_repl1$Pyruvate_repl1)
max_x_pyruvate_secretion_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
pyruvate_secretion_rate_cer_1 <- (as.numeric(pyruvate_secretion_cer_1)/Pyr_mw)*as.numeric(grate_cer_1)

# plot pyruvate secretion yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1 , y=Pyruvate_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_pyruvate_secretion_cer_1*(0.7), y=max_y_pyruvate_secretion_cer_1*(0.8), label= r_squared_pyruvate_secretion_cer_1) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_1*(0.7), y=max_y_pyruvate_secretion_cer_1*(0.9), label= pyruvate_secretion_cer_1) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_1*(0.9), y=max_y_pyruvate_secretion_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_pyruvate_secretion_cer_1*(0.9), y=max_y_pyruvate_secretion_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Pyruvate_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
pyruvate_secretion_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
pyruvate_secretion_cer_2 <- substring(pyruvate_secretion_cer_2,1,7)
r_squared_pyruvate_secretion_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_pyruvate_secretion_cer_2 <- substring(r_squared_pyruvate_secretion_cer_2,1,6)
max_y_pyruvate_secretion_cer_2 <- max(cerevisiae_repl2$Pyruvate_repl2)
max_x_pyruvate_secretion_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
pyruvate_secretion_rate_cer_2 <- (as.numeric(pyruvate_secretion_cer_2)/Pyr_mw)*as.numeric(grate_cer_2)

# plot pyruvate secretion yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2 , y=Pyruvate_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_pyruvate_secretion_cer_2*(0.7), y=max_y_pyruvate_secretion_cer_2*(0.8), label= r_squared_pyruvate_secretion_cer_2) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_2*(0.7), y=max_y_pyruvate_secretion_cer_2*(0.9), label= pyruvate_secretion_cer_2) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_2*(0.9), y=max_y_pyruvate_secretion_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_pyruvate_secretion_cer_2*(0.9), y=max_y_pyruvate_secretion_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Pyruvate_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
pyruvate_secretion_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
pyruvate_secretion_cer_3 <- substring(pyruvate_secretion_cer_3,1,7)
r_squared_pyruvate_secretion_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_pyruvate_secretion_cer_3 <- substring(r_squared_pyruvate_secretion_cer_3,1,6)
max_y_pyruvate_secretion_cer_3 <- max(cerevisiae_repl3$Pyruvate_repl3)
max_x_pyruvate_secretion_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
pyruvate_secretion_rate_cer_3 <- (as.numeric(pyruvate_secretion_cer_3)/Pyr_mw)*as.numeric(grate_cer_3)

# plot pyruvate secretion yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3 , y=Pyruvate_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_pyruvate_secretion_cer_3*(0.7), y=max_y_pyruvate_secretion_cer_3*(0.8), label= r_squared_pyruvate_secretion_cer_3) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_3*(0.7), y=max_y_pyruvate_secretion_cer_3*(0.9), label= pyruvate_secretion_cer_3) + 
  annotate("text", x=max_x_pyruvate_secretion_cer_3*(0.9), y=max_y_pyruvate_secretion_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_pyruvate_secretion_cer_3*(0.9), y=max_y_pyruvate_secretion_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# acetate secretion rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Acetate_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
acetate_secretion_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
acetate_secretion_cer_1 <- substring(acetate_secretion_cer_1,1,7)
r_squared_acetate_secretion_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_acetate_secretion_cer_1 <- substring(r_squared_acetate_secretion_cer_1,1,6)
max_y_acetate_secretion_cer_1 <- max(cerevisiae_repl1$Acetate_repl1)
max_x_acetate_secretion_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
acetate_secretion_rate_cer_1 <- (as.numeric(acetate_secretion_cer_1)/Ace_mw)*as.numeric(grate_cer_1)

# plot acetate secretion yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1 , y=Acetate_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_acetate_secretion_cer_1*(0.7), y=max_y_acetate_secretion_cer_1*(0.8), label= r_squared_acetate_secretion_cer_1) + 
  annotate("text", x=max_x_acetate_secretion_cer_1*(0.7), y=max_y_acetate_secretion_cer_1*(0.9), label= acetate_secretion_cer_1) + 
  annotate("text", x=max_x_acetate_secretion_cer_1*(0.9), y=max_y_acetate_secretion_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_acetate_secretion_cer_1*(0.9), y=max_y_acetate_secretion_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Acetate_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
acetate_secretion_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
acetate_secretion_cer_2 <- substring(acetate_secretion_cer_2,1,7)
r_squared_acetate_secretion_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_acetate_secretion_cer_2 <- substring(r_squared_acetate_secretion_cer_2,1,6)
max_y_acetate_secretion_cer_2 <- max(cerevisiae_repl2$Acetate_repl2)
max_x_acetate_secretion_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
acetate_secretion_rate_cer_2 <- (as.numeric(acetate_secretion_cer_2)/Ace_mw)*as.numeric(grate_cer_2)

# plot acetate secretion yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2, y=Acetate_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_acetate_secretion_cer_2*(0.7), y=max_y_acetate_secretion_cer_2*(0.8), label= r_squared_acetate_secretion_cer_2) + 
  annotate("text", x=max_x_acetate_secretion_cer_2*(0.7), y=max_y_acetate_secretion_cer_2*(0.9), label= acetate_secretion_cer_2) + 
  annotate("text", x=max_x_acetate_secretion_cer_2*(0.9), y=max_y_acetate_secretion_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_acetate_secretion_cer_2*(0.9), y=max_y_acetate_secretion_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Acetate_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
acetate_secretion_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
acetate_secretion_cer_3 <- substring(acetate_secretion_cer_3,1,7)
r_squared_acetate_secretion_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_acetate_secretion_cer_3 <- substring(r_squared_acetate_secretion_cer_3,1,6)
max_y_acetate_secretion_cer_3 <- max(cerevisiae_repl3$Acetate_repl3)
max_x_acetate_secretion_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
acetate_secretion_rate_cer_3 <- (as.numeric(acetate_secretion_cer_3)/Ace_mw)*as.numeric(grate_cer_3)

# plot acetate secretion yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3, y=Acetate_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_acetate_secretion_cer_3*(0.7), y=max_y_acetate_secretion_cer_3*(0.8), label= r_squared_acetate_secretion_cer_3) + 
  annotate("text", x=max_x_acetate_secretion_cer_3*(0.7), y=max_y_acetate_secretion_cer_3*(0.9), label= acetate_secretion_cer_3) + 
  annotate("text", x=max_x_acetate_secretion_cer_3*(0.9), y=max_y_acetate_secretion_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_acetate_secretion_cer_3*(0.9), y=max_y_acetate_secretion_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# succinate secretion rate
# replicate 1
y_cer_1 <- cerevisiae_repl1$Succinate_repl1
x_cer_1 <- cerevisiae_repl1$CDW_repl1
linear_model_cer_1 <- lm(y_cer_1 ~ x_cer_1)
succinate_secretion_cer_1 <- as.numeric(linear_model_cer_1$coefficients[2])
succinate_secretion_cer_1 <- substring(succinate_secretion_cer_1,1,7)
r_squared_succinate_secretion_cer_1 <- as.character(summary(linear_model_cer_1)$r.squared)
r_squared_succinate_secretion_cer_1 <- substring(r_squared_succinate_secretion_cer_1,1,6)
max_y_succinate_secretion_cer_1 <- max(cerevisiae_repl1$Succinate_repl1)
max_x_succinate_secretion_cer_1 <- max(cerevisiae_repl1$CDW_repl1)
succinate_secretion_rate_cer_1 <- (as.numeric(succinate_secretion_cer_1)/Succ_mw)*as.numeric(grate_cer_1)

# plot succinate secretion yield
ggplot(cerevisiae_repl1, aes(x=CDW_repl1, y=Succinate_repl1)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_1$coefficients[1], slope = linear_model_cer_1$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_succinate_secretion_cer_1*(0.7), y=max_y_succinate_secretion_cer_1*(0.8), label= r_squared_succinate_secretion_cer_1) + 
  annotate("text", x=max_x_succinate_secretion_cer_1*(0.7), y=max_y_succinate_secretion_cer_1*(0.9), label= succinate_secretion_cer_1) + 
  annotate("text", x=max_x_succinate_secretion_cer_1*(0.9), y=max_y_succinate_secretion_cer_1*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_succinate_secretion_cer_1*(0.9), y=max_y_succinate_secretion_cer_1*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 2
y_cer_2 <- cerevisiae_repl2$Succinate_repl2
x_cer_2 <- cerevisiae_repl2$CDW_repl2
linear_model_cer_2 <- lm(y_cer_2 ~ x_cer_2)
succinate_secretion_cer_2 <- as.numeric(linear_model_cer_2$coefficients[2])
succinate_secretion_cer_2 <- substring(succinate_secretion_cer_2,1,7)
r_squared_succinate_secretion_cer_2 <- as.character(summary(linear_model_cer_2)$r.squared)
r_squared_succinate_secretion_cer_2 <- substring(r_squared_succinate_secretion_cer_2,1,6)
max_y_succinate_secretion_cer_2 <- max(cerevisiae_repl2$Succinate_repl2)
max_x_succinate_secretion_cer_2 <- max(cerevisiae_repl2$CDW_repl2)
succinate_secretion_rate_cer_2 <- (as.numeric(succinate_secretion_cer_2)/Succ_mw)*as.numeric(grate_cer_2)

# plot succinate secretion yield
ggplot(cerevisiae_repl2, aes(x=CDW_repl2, y=Succinate_repl2)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_2$coefficients[1], slope = linear_model_cer_2$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_succinate_secretion_cer_2*(0.7), y=max_y_succinate_secretion_cer_2*(0.8), label= r_squared_succinate_secretion_cer_2) + 
  annotate("text", x=max_x_succinate_secretion_cer_2*(0.7), y=max_y_succinate_secretion_cer_2*(0.9), label= succinate_secretion_cer_2) + 
  annotate("text", x=max_x_succinate_secretion_cer_2*(0.9), y=max_y_succinate_secretion_cer_2*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_succinate_secretion_cer_2*(0.9), y=max_y_succinate_secretion_cer_2*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# replicate 3
y_cer_3 <- cerevisiae_repl3$Succinate_repl3
x_cer_3 <- cerevisiae_repl3$CDW_repl3
linear_model_cer_3 <- lm(y_cer_3 ~ x_cer_3)
succinate_secretion_cer_3 <- as.numeric(linear_model_cer_3$coefficients[2])
succinate_secretion_cer_3 <- substring(succinate_secretion_cer_3,1,7)
r_squared_succinate_secretion_cer_3 <- as.character(summary(linear_model_cer_3)$r.squared)
r_squared_succinate_secretion_cer_3 <- substring(r_squared_succinate_secretion_cer_3,1,6)
max_y_succinate_secretion_cer_3 <- max(cerevisiae_repl3$Succinate_repl3)
max_x_succinate_secretion_cer_3 <- max(cerevisiae_repl3$CDW_repl3)
succinate_secretion_rate_cer_3 <- (as.numeric(succinate_secretion_cer_3)/Succ_mw)*as.numeric(grate_cer_3)

# plot succinate secretion yield
ggplot(cerevisiae_repl3, aes(x=CDW_repl3, y=Succinate_repl3)) +
  geom_point(color="red",size=2) +
  geom_abline(intercept = linear_model_cer_3$coefficients[1], slope = linear_model_cer_3$coefficients[2]) +
  theme_bw() +
  annotate("text", x=max_x_succinate_secretion_cer_3*(0.7), y=max_y_succinate_secretion_cer_3*(0.8), label= r_squared_succinate_secretion_cer_3) + 
  annotate("text", x=max_x_succinate_secretion_cer_3*(0.7), y=max_y_succinate_secretion_cer_3*(0.9), label= succinate_secretion_cer_3) + 
  annotate("text", x=max_x_succinate_secretion_cer_3*(0.9), y=max_y_succinate_secretion_cer_3*(0.8), label= '(R2)') + 
  annotate("text", x=max_x_succinate_secretion_cer_3*(0.9), y=max_y_succinate_secretion_cer_3*(0.9), label= 'g/gDW') + 
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank())

# collect data into data frame
collected_HPLCdata_cerevisiae <- data.frame(matrix(nrow = 8,ncol = 3))
colnames(collected_HPLCdata_cerevisiae) <- c('Repl1','Repl2','Repl3')
collected_HPLCdata_cerevisiae[1,] <- data.frame((-1/as.numeric(glucose_uptake_cer_1)),(-1/as.numeric(glucose_uptake_cer_2)),(-1/as.numeric(glucose_uptake_cer_3)))
collected_HPLCdata_cerevisiae[2,] <- data.frame(glucose_uptake_rate_cer_1,glucose_uptake_rate_cer_2,glucose_uptake_rate_cer_3)
collected_HPLCdata_cerevisiae[3,] <- data.frame(glycerol_secretion_rate_cer_1,glycerol_secretion_rate_cer_2,glycerol_secretion_rate_cer_3)
collected_HPLCdata_cerevisiae[4,] <- data.frame(ethanol_secretion_rate_cer_1,ethanol_secretion_rate_cer_2,ethanol_secretion_rate_cer_3)
collected_HPLCdata_cerevisiae[5,] <- data.frame(pyruvate_secretion_rate_cer_1,pyruvate_secretion_rate_cer_2,pyruvate_secretion_rate_cer_3)
collected_HPLCdata_cerevisiae[6,] <- data.frame(acetate_secretion_rate_cer_1,acetate_secretion_rate_cer_2,acetate_secretion_rate_cer_3)
collected_HPLCdata_cerevisiae[7,] <- data.frame(succinate_secretion_rate_cer_1,succinate_secretion_rate_cer_2,succinate_secretion_rate_cer_3)
collected_HPLCdata_cerevisiae[8,] <- data.frame(as.numeric(grate_cer_1),as.numeric(grate_cer_2),as.numeric(grate_cer_3))

collected_HPLCdata_cerevisiae <- collected_HPLCdata_cerevisiae %>%
  rowwise() %>%
  mutate(meanValue = mean(c(Repl1,Repl2,Repl3))) %>%
  mutate(sd = sd(c(Repl1,Repl2,Repl3)))

collected_HPLCdata_cerevisiae <- data.frame(collected_HPLCdata_cerevisiae)
row.names(collected_HPLCdata_cerevisiae) <- c('BiomassYield','GlucoseUptakeRate','GlycerolProdRate',
                                              'EthanolProdRate','PyruvateProdRate','AcetateProdRate',
                                              'SuccinateProdRate','growthRate')
# Save results to file
#setwd('../../results/')
#write.table(collected_HPLCdata_cerevisiae,'rates_cerevisiae.txt',sep = '\t')


