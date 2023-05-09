# Introduction ----
# Title: Figure 7. Erythrocyte cycle
# Description: Quantitative simulations of the hemoglobin oxygen dissociation curve during circulation.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 10/05/2023

# Functions ----
source(file = "code/RBC_oxygen_model.R")

# Simulations ----
# Set the simulation variables.
po2 = 0:100                           # Create a vector of value for the partial pressure of oxygen.
bpg23_total <- 5 * (10 ^ (-3))        # Total 2,3-bisphosphoglycerate in the erythrocyte.
bpg23_lungs <- 0.9 * (10 ^ (-3))     # Bound 2,3-BPG concentration in the lungs.  
bpg23_muscle <- 4.65 * (10 ^ (-3))     # Bound 2,3-BPG concentration in the skeletal muscles.

# A. Lungs ----
# 1. Calculate p50 in the lungs
p50_lungs <- model_p50(dpg_rbc = bpg23_lungs, dpg_s = bpg23_total)
round(p50_lungs, 1)

# 2. Calculate the oxygen saturation in the lungs using the Hill equation.
SHbO2_lungs <- model_hill(po2 = po2, p50 = p50_lungs)

# 3. Create a dataframe using the lungs oxygen saturation and partial pressure of oxygen columns.
lung_data <- data.frame(po2, SHbO2_lungs)

# B. Muscles ----
# 1. Calculate p50 in the skeletal muscle.
p50_muscle <- model_p50(dpg_rbc = bpg23_muscle, dpg_s = bpg23_total)
round(p50_muscle, 1)

# 2. Calculate the oxygen saturation in the muscle with the Hill equation.
SHbO2_muscle <- model_hill(po2 = po2, p50 = p50_muscle)

# 3. Create a dataframe using the skeletal muscle oxygen saturation and partial pressure of oxygen columns.
muscle_data <- data.frame(po2, SHbO2_muscle)

# Plots ----
plot_hill(p50 = p50_lungs, 
          p50_std = p50_muscle, 
          add.std.p50 = T, 
          add.new.p50 = T, 
          add.arrow = T,
          add.text = T)

# Export ----
# Export the simulated data to ".csv" files.
write.csv2(lung_data, file = "data/lungs_data.csv", 
          row.names = FALSE)

write.csv2(muscle_data,"data/muscle_data.csv",
           row.names = FALSE)
