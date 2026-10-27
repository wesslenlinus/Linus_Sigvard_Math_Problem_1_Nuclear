import math as math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

NA = 6.022 * 10 ** 23 # Avogadro's number, pcs/mol

# Technical specifications for Ringhals 1, all are taken from the official specifications from Vattenfall
reactor_power_wh = 2540e6         
core_assemblies = 648
assembly_mass_UO2_g = 193e3    # 193 kg UO2 per assembly (public data from Vattenfall)
core_height_m = 3.7                 
core_diameter_m = 5.95             
core_mass = core_assemblies * assembly_mass_UO2_g
reactor_type = "BWR"  

# fuel specifications
fuel_material = "UO2"  
Ma_UO2 = 270  # g/mol UO2
rho_UO2 = 10.97  # g/cm^3             
enrichment_percent = 3/100  # Assumption based on course material
N_U = (rho_UO2 * NA) / Ma_UO2
N_U235 = enrichment_percent * N_U
N_U238 = (1 - enrichment_percent) * N_U 

# moderator specifications
moderator = "H2O"                  
moderator_temp_c = (286.0 + 273)    
void_fraction = 0.40                
rho_steam = 36.53e-3  # g/cm^3 at T = 286 C and P = 7 MPa
rho_water = 739.7e-3  # g/cm^3
rho_mod = void_fraction * rho_steam + (1 - void_fraction) * rho_water
Ma_mod = 18  # g/mol for H2O       
mod_fuel_volume_ratio = 3        # (Vm/Vf) assumed given lectures in course
sigma_I_abs = 20.0 * 10**-24    # barn converted to cm^2
xi_Sigma_s = 0.8                # cm^-1 (assumed value for H2O)

# cross sections from nndc.bnl.gov
sigma_f_235 = 596 * 10 ** -24  # barn
sigma_a_235 = 97.6 * 10 ** -24 + sigma_f_235  #(n,gamma) + (n, total fission)  # barn
sigma_f_238 = 1.76 * 10 ** -5 * 10 ** -24
sigma_a_238 = 2.78 * 10 ** -24 + sigma_f_238

# Fissile and fertile cores at the start of the cycle
energy_per_fission = 3.2e-11     # Joules / Fission
operating_days = 365 * 1.5
energy_per_cycle = reactor_power_wh * 60 * 60 * 24 * operating_days  # Assumption: 1 cycle = 18 months
N_U235_start = enrichment_percent * core_mass * NA / Ma_UO2  

# Fissile and fertile cores at the end of the cycle
fissions_per_cycle = energy_per_cycle / energy_per_fission 
fissile_gone = fissions_per_cycle * 0.95 * (sigma_a_235 / sigma_f_235)
N_U235_end = N_U235_start - fissile_gone
fuel_left_end = core_mass * NA / Ma_UO2
enrichment_percent_end = N_U235_end / fuel_left_end

# calculation of thermal fission factor n
def calc_n(nu_bar, e, sigmaF235, sigmaC235, sigmaC238):
    N235 = e * N_U
    N238 = (1 - e) * N_U
    n = nu_bar * (e * sigmaF235) / (e * (sigmaF235 + sigmaC235) + (1 - e) * sigmaC238)
    return n

def calc_p(e_frac, N_uran, sigma_I, xi_sigma_s, vm_vf_ratio):
    exponent = -(1.0 - e_frac) * (N_uran * sigma_I / xi_sigma_s) * (1.0 / vm_vf_ratio)
    p = (math.e ** (exponent))
    return p

# Calculating macroscopic absorption cross section for moderator (H2O)
sigma_a_H = 0.332e-24    #Reactor Physics KSU
sigma_a_H2O = 2 * sigma_a_H     #Presuming absorption for Oxygen = 0 (Reactor Physics KSU)
Sigma_a_H2O = (rho_mod * NA) / Ma_mod * sigma_a_H2O

# Constants for Reactor poisoning from 135Xe
sigma_a_Xe = 2.6e-18 
yield_Xe = 1.016e-3     
yield_I = 6.437e-2      
T_half_Xe =  33120
Lambda_Xe = math.log(2) / T_half_Xe

# Next step is calculating the safety rods physical specifications. In Ringhals the safety rods are B4C.
rho_B4C = 2.52  #g/cm^3
Ma_B4C = 55.26  #g/mol
number_of_rods = 157    #From Ringhals specifications
blade_width_rod = 15   #cm Assumed value based on course material
thickness_rod = 0.5  #cm Assumed value based on course material
N_B4C = rho_B4C * NA / Ma_B4C
sigma_a_B10 = 3840e-24
B10_enrichment = 0.199  # Assuming natural boron is ~19.9% B-10

# 4 Boron atoms per B4C molecule, scaled by B-10 enrichment
Sigma_a_bor = 4 * N_B4C * sigma_a_B10 * B10_enrichment

def calc_mod_and_rod_volumes(rod_insertion, scale_factor):
    length_rod = core_height_m * 100 * rod_insertion 
    V_fuel = core_mass / rho_UO2    
    
    # Corrected area: (width*thickness)*2 - (thickness^2) to avoid double counting the center cross
    cross_area = (2 * blade_width_rod * thickness_rod) - (thickness_rod**2)
    V_tot_rods = number_of_rods * (length_rod * cross_area)
    
    # Calibrated by making sure k_eff = 1 at start. 
    V_rel_rods = V_tot_rods / V_fuel * scale_factor

    # We now enforce the assumed volume ratio for the moderator (3) instead of geometric calculation
    V_rel_mod = mod_fuel_volume_ratio 
    
    return V_rel_rods, V_rel_mod

# ε (epsilon) - Fast fission factor
epsilon = 1.06

# P - Total non-leakage factor
P_leakage = 0.97

def calc_k_eff(n, E, p, f, P):
    return n * E * p * f * P

# Making sure k = 1 at start and end by varying rod insertion. 
def evaluate_reactor_state(rod_insertion, current_enrichment, current_void=0.40):
    
    # 1. Calculate moderator density with void as a variable. 
    current_rho_mod = current_void * rho_steam + (1 - current_void) * rho_water
    density_ratio = current_rho_mod / (0.40 * rho_steam + (1 - 0.40) * rho_water)
    
    # Rescale parameters for absorption and moderating based on density
    current_xi_Sigma_s = xi_Sigma_s * density_ratio
    current_Sigma_a_H2O = Sigma_a_H2O * density_ratio

    # Calculate n and p with current enrichment
    n_val = calc_n(2.43, current_enrichment, sigma_f_235, sigma_a_235-sigma_f_235, sigma_a_238-sigma_f_238)
    p_val = calc_p(current_enrichment, N_U, sigma_I_abs, current_xi_Sigma_s, mod_fuel_volume_ratio)
    
    # Find the fuels microscopic absorption cross section as enrichment varies
    N_235_current = current_enrichment * N_U
    N_238_current = (1 - current_enrichment) * N_U
    sig_a_fuel = sigma_a_235 * N_235_current + sigma_a_238 * N_238_current
    
    # Dynamic Xenon calculation based on current enrichment and flux
    sig_f_fuel_current = sigma_f_235 * N_235_current
    N_U235_core_current = current_enrichment * core_mass * NA / Ma_UO2
    current_neutron_flux = reactor_power_wh / (sigma_f_235 * energy_per_fission * N_U235_core_current)
    
    N_Xe_current = (yield_I + yield_Xe) * sig_f_fuel_current * current_neutron_flux / (Lambda_Xe + sigma_a_Xe * current_neutron_flux)
    sig_a_Xe_current = N_Xe_current * sigma_a_Xe
    
    # Calculate volumes and (f) as the rod insertion changes
    V_rel_rods, V_rel_mod = calc_mod_and_rod_volumes(rod_insertion, scale_factor=0.15)
    f_val = sig_a_fuel / (sig_a_fuel + current_Sigma_a_H2O * V_rel_mod + sig_a_Xe_current + Sigma_a_bor * V_rel_rods)
    
    k_eff = calc_k_eff(n_val, epsilon, p_val, f_val, P_leakage)
    
    return k_eff, n_val, epsilon, p_val, f_val, P_leakage

# Find rod insertion for k = 1
def target_critical(rod_insertion, current_enrichment):
    return evaluate_reactor_state(rod_insertion, current_enrichment)[0] - 1.0

rod_critical_start = brentq(target_critical, 0.0, 0.99, args=(enrichment_percent,))
rod_critical_end = brentq(target_critical, 0.0, 0.99, args=(enrichment_percent_end,))

# --- EXTRACT AND PRINT FACTORS FOR START AND END ---
k_start, n_start, eps_start, p_start, f_start, P_start = evaluate_reactor_state(rod_critical_start, enrichment_percent)
k_end, n_end, eps_end, p_end, f_end, P_end = evaluate_reactor_state(rod_critical_end, enrichment_percent_end)

print("FEMFAKTORFORMELN & REAKTORSTATUS")
print(f"{'Parameter':<10}  {'BOC (Start)':<12} {'EOC (Slut)':<12}")
print(f"{'k_eff':<10} {k_start:<12.5f} {k_end:<12.5f}")
print(f"{'eta (n)':<10} {n_start:<12.5f} {n_end:<12.5f}")
print(f"{'epsilon':<10} {eps_start:<12.5f} {eps_end:<12.5f}")
print(f"{'p':<10} {p_start:<12.5f} {p_end:<12.5f}")
print(f"{'f':<10} {f_start:<12.5f} {f_end:<12.5f}")
print(f"{'P':<10} {P_start:<12.5f} {P_end:<12.5f}")
print(f"Styrstavar {rod_critical_start*100:<11.1f}% {rod_critical_end*100:<11.1f}%")

# Plotting void coefficients
void_fractions = np.linspace(0.0, 0.8, 50)
k_eff_start_list = []
k_eff_end_list = []

for v in void_fractions:
    k_start_val = evaluate_reactor_state(rod_critical_start, enrichment_percent, current_void=v)[0]
    k_end_val = evaluate_reactor_state(rod_critical_end, enrichment_percent_end, current_void=v)[0]
    k_eff_start_list.append(k_start_val)
    k_eff_end_list.append(k_end_val)

plt.figure(figsize=(10, 6))
plt.plot(void_fractions * 100, k_eff_start_list, label="Start of Cycle (BOC)", color="blue", linewidth=2)
plt.plot(void_fractions * 100, k_eff_end_list, label="End of Cycle (EOC)", color="red", linewidth=2, linestyle='--')
plt.axvline(40, color='gray', linestyle=':', label='Operational Void (40%)')
plt.axhline(1.0, color='black', linestyle='-', linewidth=1, label='Criticality (k = 1)')

plt.title("Void Coefficient: $k_{eff}$ vs Void Fraction", fontsize=14)
plt.xlabel("Void Fraction (%)", fontsize=12)
plt.ylabel("Effective Multiplication Factor ($k_{eff}$)", fontsize=12)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()
