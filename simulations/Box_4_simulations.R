# Introduction ----
# Title: Modeling erythrocyte oxygen delivery: Box 4
# Purpose: Numerical simulations of total oxygen delivery capacity of blood
#
# Author: Panagiotis Chatzinikolaou
# Contact: chatzinpn@phed-sr.auth.gr
#
# Latest update: 27/10/2023
#
# Source: Oxygen transport: A Redox O2dyssey

compute_oxygen_delivery(co = 5, total.oxygen = compute_oxygen_capacity_2(hb = 150, bind = 1.34, sat = 0.97))

compute_oxygen_delivery(co = 5, total.oxygen = compute_oxygen_capacity_2(hb = 160, bind = 1.34, sat = 0.97))

compute_oxygen_delivery(co = 30, total.oxygen = compute_oxygen_capacity_2(hb = 150, bind = 1.34, sat = 0.97))

compute_oxygen_delivery(co = 30, total.oxygen = compute_oxygen_capacity_2(hb = 160, bind = 1.34, sat = 0.97))
