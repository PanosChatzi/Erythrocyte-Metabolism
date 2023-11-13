# Introduction ----
# Title: Erythrocyte capillary transit time
#
# Description: Mathematical equations to estimate erythrocyte capillary transit time.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 04/11/2023

# Equations ----

# Function to calculate the individual capillary cross-sectional area from its diameter.
compute_cap_area <- function(cap_diameter) { # mean capillary diameter (μm)
  
  cap_radius <- cap_diameter / 2             # individual capillary radius (μm)
  
  cap_area <- pi * (cap_radius^2)            # individual capillary CSA (μm^2)
  
  return(cap_area)
} 

# Define the function to calculate capillary transit time.
compute_transit_time <- function(cap_area,        # individual capillary CSA (μm^2)
                                 cap_density,     # capillary density (capillaries / mm^2)
                                 blood_flow_rate, # remember to convert to mL/s per g
                                 tortuosity_constant = 1.2) {
  
  # Step 1. Calculate total capillary area
  total_cap_area <- cap_area * cap_density
  
  # Step 2. Calculate the capillary fractional area per fiber area (1 mm^2 or 1.000.000 μm^2)
  fiber_area <- 1000000
  fractional_area_per_fiber <- total_cap_area / fiber_area
  
  # Step 3. Adjust for tortuosity and branching using a constant of 1.2 
  # (1.0 = perfect anisotropy & 2.0 = random orientation).
  adjusted_area_per_fiber <- fractional_area_per_fiber * tortuosity_constant
  
  # Step 4. Divide by muscle blood_flow_rate (mL/s per g).
  transit_time <- adjusted_area_per_fiber / blood_flow_rate
  
  return(transit_time)
}

# Assumptions of the function based on Richardson et al (1994).
# The 4th step assumes that leg cardiac output (Q) is homogeneously distributed within and between capillaries in the exercising muscle, 
# utilizes the relationship between capillaries per fiber volume and the muscle volume normalized leg Q to determine red cell transit time.

# Calculations ----
# Firstly, we calculate the individual capillary volume
single_cap_vol <- compute_cap_area(cap_diameter = 6)

# Secondly, we derive the leg blood flow value from Table 1 (Richardson et al, 1993)
# The leg blood flow during 100% of work is 9.10 L/min.

# To calculate the mass specific blood flow (L/min per kg) we divide it with the 
# average quadriceps femoris mass weight (2.36 kg) => 9.10/2.36 = 3.85 L/min/kg

# Thirdly, we have to convert the 3.85 L/min per kg to mL/s per g
# 1 L = 1000 mL
# 1 min = 60 s
# 1 kg = 1000 g

# 3.85 L/min per kg = 3.85 * 1000 mL/min per kg
# 3.85 L/min per kg = 3850 mL/min per kg
# 3.85 L/min per kg = 3850 mL/(60 s)/(1000 g)

blood_flow <- 3850 / 60 / 1000 # 3850 mL/(60 s)/(1000 g) = 0.06416667 mL/s per g

# Finally, we compute the capillary transit time using the custom function (= 0.15863 s)
compute_transit_time(cap_area = single_cap_vol, cap_density = 300, blood_flow_rate = blood_flow)

# Using the values employed in Richardson et al. (1994) to calculate a transit time of ≈0.11 s
compute_transit_time(cap_area = 28.3, cap_density = 200, blood_flow_rate = 0.0641667)

# References
# The calculations were derived from:
#
# 1. Richardson RS et al. Red blood cell transit time in man: theoretical effects of capillary density. 
#    Adv Exp Med Biol. 1994;361:521-32. doi: 10.1007/978-1-4615-1875-4_91. PMID: 7597979.
#
# 2. Richardson RS et al. High muscle blood flow in man: is maximal O2 extraction compromised? 
#    J Appl Physiol (1985). 1993 Oct;75(4):1911-6. doi: 10.1152/jappl.1993.75.4.1911. PMID: 8282650.