# Introduction ----
# Title: Erythrocyte Metabolism 
# Description: Mathematical models of the hemoglobin oxygen and carbon dioxide dissociation curves.
#
# Author: Panagiotis N. Chatzinikolaou
# Affiliation: PhD candidate at Aristotle University of Thessaloniki, Greece
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 09/05/2023

# Oxygen curve models ----
# Custom functions based on the equations and data of Dash et al, 2016.

# Calculate the Hill coefficient nH using equation 11
compute_nHill <- function(po2 = 100) { # Partial pressure of oxygen in standard conditions
  
  # Assign the parameter values based on Table 1
  alpha <- 2.82
  beta <- 1.20
  gamma <- 29.25

  nH <-  alpha - beta * (10 ^ (- (po2 / gamma)))
  nH
}

# Hill model based on equations 6 and 11 from Dash et al, 2016.
model_hill <- function(po2 = 100,     # partial oxygen pressure in standard conditions.
                       p50 = 26.8,    # p50 in standard conditions.
                       nhill.var = TRUE) {
  
  # If nHill coefficient is dynamic, then calculate it using equation 11
  # of Dash, Korman & Bassingthwaighte 2016
  if(nhill.var == TRUE) {
    
    nH <-  compute_nHill(po2 = po2)
  
  } else {
    
    # If nHill.var is FALSE, then assign the standard value.
    nH <- 2.7
  }
  
  # Calculate the oxygen saturation of hemoglobin
  x <- ((po2 / p50) ^ nH) / (1 + (po2 / p50) ^ nH)
  
  return(x)
}

# Function to calculate p50 (mmHg) based on pH, pCO2, 2,3-BPG and temperature.
# Equations 9a-9d and equation 10. Data derived from Table 1.
model_p50 <- function(pH_rbc = 7.24,                # pH standard conditions
                      pco2 = 40,                    # pco2 standard conditions
                      dpg_rbc = 4.65 * (10 ^ (-3)), # 2,3-BPG standard conditions
                      temp = 37,                    # Temperature standard conditions
                      p50_s = 26.8,                 # p50 in standard conditions
                      pH_s = 7.24,                  # pH standard conditions
                      pco2_s = 40,                  # pco2 standard conditions
                      dpg_s = 4.65 * (10 ^ (-3)),   # 2,3-BPG standard conditions
                      temp_s = 37) {                  
  
  p50dpH <-  p50_s - 25.535 * (pH_rbc - pH_s) + 10.646 * ((pH_rbc - pH_s) ^ 2) -
    1.764 * ((pH_rbc - pH_s) ^ 3)
  
  p50dco2 <-  p50_s + 0.1273 * (pco2 - pco2_s) + 1.083 * (10 ^ (- 4)) * ((pco2 - pco2_s) ^ 2)
  
  p50dBPG <-  p50_s + 795.63 * (dpg_rbc - dpg_s) - 19660.89 * ((dpg_rbc - dpg_s) ^ 2)
  
  p50dtemp <-  p50_s + 1.435 * (temp - temp_s) + (4.163 * (10 ^ (- 2))) * ((temp - temp_s) ^ 2) +
    (6.86 * (10 ^ (- 4))) * ((temp - temp_s) ^ 3)
  
  # Calculate the p50 based on equation 10
  p50 <- p50_s * (p50dpH / p50_s) * (p50dco2 / p50_s) * (p50dBPG / p50_s) * (p50dtemp / p50_s)
  
  return(p50)
}            # Temperature in standard conditions

# Erythrocyte pH function Based on Siggaard-Andersen and colleagues (1971)
rbc_ph <- function (ph = 7.4) { # Plasma pH in standard conditions
  
  y <- 0.795 * ph + 1.357
  
  return(y)
}

# Function to calculate free oxygen in the water space of erythrocytes in M/mmHg (Equation 2a).
oxy_bl <- function(po2 = 100,      # Partial oxygen pressure in standard conditions
                   temp = 37,      # Temperature in standard conditions
                   w_pl = 0.94) {  # Fractional water space of plasma
  
  ao2 <- (1.37 - 0.0137 * (temp - 37) + 0.00058 * 
            ((temp - 37) ^ 2)) * ((10 ^ (-6)) / w_pl)
  
  x <- ao2 * po2
  
  return(x)
}

# Function to calculate free carbon dioxide in the water space of erythrocytes in M/mmHg (Equation 2b).
co2_bl <- function(pco2 = 40,       # Partial carbon dioxide pressure in standard conditions
                   temp = 37,       # Temperature in standard conditions
                   w_pl = 0.94) {   # Fractional water space of plasma
  
  aco2 <- (3.07 - 0.057 * (temp - 37) + 0.002 * 
             (temp - 37) ^ 2) * ((10 ^ (-5)) / w_pl)
  
  x <- aco2 * pco2
  
  return(x)
}

# Function to calculate F1, F2, F3, F4 (Equations 4a-d)
compute_F <- function (ph = 7.4) { # Plasma pH in standard conditions.
  
  # First, calculate erythrocyte pH.
  pHrbc <- 0.795 * ph + 1.357  

  # Second, calculate hydrogen concentration in erythrocytes. 
  hydrogen_con <- 10 ^ (- pHrbc)

  # Third, assign standard parameters based on Table 1.
  K2_prime2 <- 10 ^ -6
  K3_prime2 <- 10 ^ -6
  K5_prime2 <- 2.4 * (10 ^ -8)
  K6_prime2 <- 1.2 * (10 ^ -8)
  
  # Fourth, calculate F1, F2, F3 and F4.
  F1 <- 1 + (K2_prime2 / hydrogen_con)

  F2 <- 1 + (K3_prime2 / hydrogen_con)

  F3 <- 1 + (hydrogen_con / K5_prime2)

  F4 <- 1 + (hydrogen_con / K6_prime2)
  
  # Fifth, combine all F values into a vector.
  Fvec <- c(F1, F2, F3, F4)
  
  return(Fvec)
}

# Model of HbO2 saturation (SHbO2; equation 1a) based on Dash et al, 2016.
model_dash <- function(dpg_rbc = 4.65 * (10 ^ (-3)), # 2,3-BPG standard conditions.
                       po2 = 100,       # Partial pressure of oxygen in standard conditions.
                       pco2 = 40,       # Partial pressure of carbon dioxide in standard conditions.
                       ph = 7.4,        # Plasma pH in standard conditions.
                       temp = 37,       # Temperature in standard conditions.
                       w_pl = 0.94,     # Fractional water space of plasma.
                       p50.var = TRUE,  # Evaluate if p50 is in standard conditions.
                       k4.var = TRUE) { # Evaluate if K4' in standard conditions.
  
  # This function is organized in 8 steps to calculate  oxygen saturation:
  
  # Step 1. Set the parameters based on Table 1.
  K2_prime1 <- 21.5 # or 23.65
  K3_prime1 <- 11.3 # or 14.7
  
  # Step 2. Calculate the free oxygen and carbon dioxide with equations 2a & 2b.
  free_oxy_bl <- oxy_bl(po2 = po2, temp = temp)
  free_co2_bl <- co2_bl(pco2 = pco2, temp = temp)
  
  # Step 3. Calculate F1-F4 using equations 4a-d.
  F1 <- compute_F(ph = ph)[1]
  F2 <- compute_F(ph = ph)[2]
  F3 <- compute_F(ph = ph)[3]
  F4 <- compute_F(ph = ph)[4]
  
  # Step 4. Calculate p50 using custom functions of equations 9a-d and 10.
  # Note: that way, we make SHbO2 dependent on all physiological variables,
  # namely 2,3-BPG, temperature, pH and PCO2.
  if(p50.var == TRUE) {
    
    ph_rbc <- 0.795 * ph + 1.357
    
    p50 <- model_p50(pH_rbc = ph_rbc, temp = temp,
                     pco2 = pco2, dpg_rbc = dpg_rbc)
    
  } else {
    
    # If p50 is in standard conditions.
    p50 <- 26.8
  }
  
  # Step 5. Evaluate if k4.var is TRUE.
  if(k4.var == TRUE) {
    #  1. If TRUE, the calculate K4_prime1 based on equation 7.
    
    #  2. Calculate nH using equation 11.
    nH <-  compute_nHill(po2 = po2)
    
    #  3. Calculate ao2 using equation 2a.
    ao2 <- (1.37 - 0.0137 * (temp - 37) + 0.00058 * 
              ((temp - 37) ^ 2)) * ((10 ^ (-6)) / w_pl)
    
    #  4. Calculate aco2 using equation 2b.
    aco2 <- (3.07 - 0.057 * (temp - 37) + 0.002 * 
               (temp - 37) ^ 2) * ((10 ^ (-5)) / w_pl)
    
    #  5. Now calculate K4prime.
    K4_prime1 <- (((ao2 * po2) ^ (nH - 1)) * ((K2_prime1 * aco2 * pco2 * F1) + F3)) /
      (((ao2 * p50) ^ nH) * ((K3_prime1 * aco2 * pco2 * F2) + F4))
    
  } else {
    
    # If false then assign the standard value from Table 1 (Dash et al 2016).
    K4_prime1 <- 2.03 * (10 ^ 5)
  }
  
  # Step 6. Calculate KHbO2 using equation 3a.
  khbo2 <- (K4_prime1 * ((K3_prime1 * free_co2_bl * F2) + F4)) / ((K2_prime1 * free_co2_bl * F1) + F3)
  
  # Step 7. Calculate SHbO2 using equation 1a.
  sat <-  (khbo2 * free_oxy_bl) / (1 + (khbo2 * free_oxy_bl))
  
  # Step 8. Combine the saturation, K4_prime, and khbo2 values into a list.
  out <- list(sat, K4_prime1, khbo2)
  
  names(out) <- c("SHbo2", "K4prime", "KHbo2")
  
  return(out)
}

# Oxygen delivery models ----
# Function to calculate the oxygen delivery capacity. 
compute_total_oxygen <- function(bind = 1.39,     # Binding capacity can vary between 1.30, 1.34 and 1.39.
                                 hb = 15,         # Hemoglobin concentration in grams per L.
                                 sat = 0.972,     # Oxygen saturation.
                                 po2 = 100) {     # Partial pressure of oxygen in arterial blood.
  
  # Equation from Dunn 2016
  oxygen <- (bind * hb * sat) + (0.0225 * (po2 / 7.50062))
  
  # Equation from Deranged Physiology website.
  #oxygen <- (bind * 150 * 0.972) + (0.03 * po2)
  
  return(oxygen)
}

# Based on the oxygen equation in Table 1 in Dash et al, 2016. 
compute_oxygen_dash <- function(po2 = 100,               # Partial pressure of oxygen.
                                hct = 0.45,              # Hematocrit.
                                sat = 0.972,             # Oxygen saturation of hemoglobin.
                                blood_hb = 150,          # Hemoglobin concentration in blood (mg / L).
                                temp = 37,               # Temperature in standard conditions.
                                water_plasma = 0.94,     # Fractional water space of plasma.
                                water_rbc = 0.65) {      # Fractional water space of erythrocytes.

  # 1. Fractional water space of blood.
  blood_w <- (1 - hct) * water_plasma + hct * water_rbc
  
  # 2. Solubility of oxygen in water at temperature = temp.
  ao2 <- (1.37 - 0.0137 * (temp - 37) + 0.00058 * 
            ((temp - 37) ^ 2)) * ((10 ^ (-6)) / water_plasma)
  
  # 3. Hemoglobin concentration in erythrocytes.
  hb <- blood_hb / 64458
  erythro_hb <- hb / hct
  
  # 4. Total oxygen content in blood.
  total_oxygen <- blood_w * ao2 * po2 + 4 * hct * erythro_hb * sat
  
  # 5. Convert oxygen mol/l (M) to mL oxygen per 100 mL blood.
  converting_factor <- 22.256 # (L of gas/mol)
  converted.o2 <- total_oxygen * converting_factor * 100
  
  return(converted.o2)
}

# Oxygen delivery capacity in blood.
oxygen.delivery <- function(co = 5,       # Cardiac out in L/min.
                            total.oxygen) {
  
  do2 <- co * total.oxygen
  
  return(do2)
}

# Plot functions ----
# Function to plot the oxygen dissociation curve.
plot_hill <- function(x = 0:100,            # Create a vector of values for PO2.
                      p50 = 26.8,           # p50 dynamic.
                      p50_std = 26.8,       # p50 in standard conditions.
                      add.std.p50 = FALSE,
                      add.new.p50 = FALSE,
                      add.text = FALSE,
                      add.arrow = FALSE,
                      ...) { # Optional further arguments passed on the plot.
  
  # Calculate SHbO2 based on the defined range of PO2 values.
  SHbO2 <- model_hill(po2 = x, p50 = p50) * 100
  
  # Generate the plot with optional arguments.
  plot(SHbO2,
       type = "l", lty = 1, lwd = 2,
       main= "Oxygen dissociation curve",
       xlab = "Partial pressure of oxygen (mmHg)",
       ylab = "Oxygen saturation (%)",
       xlim = c(0, 100), ylim = c(0, 100),
       las = 1, xaxs = "i", yaxs = "i",
       font.lab = 2, font.axis = 2,
       ...)
  
  # If add.std.50 is TRUE add a vertical line for the standard p50 (26.8).
  if(add.std.p50 == TRUE) {
  abline(v = p50_std, col = "black", lty = 2, lwd = 1)
    }
  
  # If add.new.50 is TRUE add a vertical line for the calculated p50.
  if(add.new.p50 == TRUE) {
  abline(v = p50, col = "red", lty = 2, lwd = 1)
    }

  # If add.text is TRUE add text showing the calculated p50 value.
  if(add.text == TRUE) {
  text(SHbO2[p50], SHbO2[p50] + 3, paste0("p50= ", round(p50, 2), " mmHg"))
    }

  # If add.arrow is TRUE add arrow to the direction of p50.
  if(add.arrow == TRUE) {
  
  # Change the direction of the arrow if p50 shifted to the left or right.
  if(p50 > p50_std) {
    arrows(x0 = p50_std - 2, y0 = 85,
           x1 = p50, y1 = 85,
           length = .075, lwd = 1.75,
           col = "red")
  } else {
    arrows(x0 = p50_std +2, y0 = 85,
           x1 = p50, y1 = 85,
           length = .075, lwd = 1.75,
           col = "red")
    }
  }
}

# References ---- 
# 1) Dash & Bassingthwaighte. 2010. (PMID: 20162361)
# 2) Dash et al. 2016. (PMID: 26298270)
# 3) Dunn. 2016. https://academic.oup.com/bjaed/article/16/10/341/2288629
# 4) Deranged Physiology website https://bit.ly/42kDXbm