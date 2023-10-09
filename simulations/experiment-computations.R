# Introduction ----
# Title: Oxygen transport computation
#
# Description: Mathematical simulations of oxygen transport using the concentrations of
# 2,3-bisphosphoglycerate before and after exercise during each experimental condition.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 09/10/2023

# Functions ----
source(file = "code/RBC_oxygen_model.R")

# Data ----
# Set the standard parameters
po2 <- 0:100  # Create a vector of values for the partial pressure of oxygen.

# Assign the experimental mean 2,3-BPG values at rest, pre- and post-exercise
# Control condition
bpg23_con_rest <- 5.05 * (10 ^ (-3))
bpg23_con_pre <- 5.11 * (10 ^ (-3))
bpg23_con_post <- 5.81 * (10 ^ (-3))

hematocrit_con <- 0.429
hemoglobin_con <- 149.15

# Oxidative stress condition
bpg23_exp_rest <- 5.08 * (10 ^ (-3))
bpg23_exp_pre <- 5.15 * (10 ^ (-3))
bpg23_exp_post <- 6.28 * (10 ^ (-3))

hematocrit_exp <- 0.434
hemoglobin_exp <- 149.95

bpg23_S <- (bpg23_con_rest + bpg23_exp_rest)/2

# Simulations ----
# Step 1. Calculate p50 values
p50_con_post <- model_p50(dpg_rbc = bpg23_con_post, dpg_s = bpg23_S)

p50_exp_post <- model_p50(dpg_rbc = bpg23_exp_post, dpg_s = bpg23_S)

# Step 2. Calculate the oxygen saturation in the muscle with the Hill equation.
SHbO2_con_post <- model_hill(po2 = 100, p50 = p50_con_post)
SHbO2_con_range <- model_hill(po2 = po2, p50 = p50_con_post)

SHbO2_exp_post <- model_hill(po2 = 100, p50 = p50_exp_post)
SHbO2_exp_range <- model_hill(po2 = po2, p50 = p50_exp_post)

# Step 3. Calculate total oxygen transport in blood for each condition
oxygen_con <- compute_oxygen_dash(sat = SHbO2_con_post, blood_hb = hemoglobin_con, hct = hematocrit_con)

oxygen_exp <- compute_oxygen_dash(sat = SHbO2_exp_post, blood_hb = hemoglobin_con, hct = hematocrit_con)

# Step 4. Create a dataframe using the SHbO2 and PO2 columns.
sim_oxy_data <- data.frame(po2, SHbO2_con_range, SHbO2_exp_range)

# Step 5. Plot the oxygen dissociation curve of the simulated data
plot_hill(p50_pre = p50_con_post, 
          p50_post = p50_exp_post,
          add.second.curve = T,
          add.std.p50 = T, 
          add.new.p50 = T,
          add.arrow = T)

# Step 6. Export the simulated data to ".csv" files.
write.csv2(sim_oxy_data, file = "data/sim_data.csv", 
           row.names = FALSE)