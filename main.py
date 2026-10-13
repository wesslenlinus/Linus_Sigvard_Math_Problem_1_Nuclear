import math as math
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

NA = 6.022 * 10 ** 23  # Avogadro's number, pcs/mol

#Technical specifications for Ringhals 1
reactor_power_wh = 2540e6         
core_assemblies = 648
assemblies_mass = 193e3                
core_height_m = 3.7                 
core_diameter_m = 5.95             
core_mass = core_assemblies * assemblies_mass
reactor_type = "BWR"  


#fuel specifications
fuel_material = "UO2"  
Ma_UO2 = 270  # g/mol UO2
rho_UO2 = 10.97  # g/cm^3             
enrichment_percent = 3/100  #Assumption based on course material
N_U = (rho_UO2 * NA) / Ma_UO2
N_U235 = enrichment_percent * N_U
N_U238 = (1 - enrichment_percent) * N_U 


#moderator specifications
moderator = "H2O"                  
moderator_temp_c = (286.0 + 273)    
void_fraction = 0.40                
rho_steam = 36.53e-3  # g/cm^3 at T = 286 C and P = 7 MPa from
rho_water = 739.7e-3  # g/cm^3
rho_mod = void_fraction * rho_steam + (1 - void_fraction) * rho_water
Ma_mod = 18  # g/mol for H2O       
mod_fuel_volume_ratio = 3        #(Vm/Vf) assumed given lectures in course
sigma_I_abs = 20.0 * 10**-24    # barn converted to cm^2
xi_Sigma_s = 0.8                # cm^-1 (assumed value for H2O)


# cross sections from nndc.bnl.gov
sigma_f_235 = 596 * 10 ** -24  # barn
sigma_a_235 = 97.6 * 10 ** -24 + sigma_f_235  #(n,gamma) + (n, total fission)  # barn
sigma_f_238 = 1.76 * 10 ** -5 * 10 ** -24
sigma_a_238 = 2.78 * 10 ** -24 + sigma_f_238


#Fissile and fertile cores at the start of the cycle
energy_per_fission = 3.2e-11     # Joules / Fission
operating_days = 365 * 1.5
energy_per_cycle = reactor_power_wh * 60 * 60 * 24 * operating_days  # Assumption: 1 cycle = 18 year
N_U235_start = enrichment_percent * core_mass * NA / Ma_UO2  #Number of fissile cores at the start of the cycle. 
                                                             #Assumption: Only fissile cores are U-235
N_U238_start = (1 - enrichment_percent) * core_mass * NA / Ma_UO2                                                             


#Fissile and fertile cores at the end of the cycle
fissions_per_cycle = energy_per_cycle / energy_per_fission 
nuclei_fissioned = fissions_per_cycle * 0.95 * (sigma_a_235 / sigma_f_235)
N_U235_end = N_U235_start - fissions_per_cycle
fuel_left_end = core_mass * NA / Ma_UO2
enrichment_percent_end = N_U235_end / fuel_left_end


#calculation of thermal fission factor n
def calc_n(nu_bar, e, sigmaF235, sigmaC235, sigmaC238):

    # Number of U235 nuclides per cm^3 fuel
    N235 = e * N_U
    # Number of U238 nuclides per cm^3 fuel
    N238 = (1 - e) * N_U

    n = nu_bar * (e * sigmaF235) / (e * (sigmaF235 + sigmaC235) + (1 - e) * sigmaC238)

    return(n)


def calc_p(e_frac, N_uran, sigma_I, xi_sigma_s, vm_vf_ratio):
    exponent = -(1.0 - e_frac) * (N_uran * sigma_I / xi_sigma_s) * (1.0 / vm_vf_ratio)
    p = (math.e ** (exponent))
    return p


thermal_fission_factor_start = calc_n(2.43, enrichment_percent, sigma_f_235, sigma_a_235-sigma_f_235, sigma_a_238-sigma_f_238)
thermal_fission_factor_end = calc_n(2.43, enrichment_percent_end, sigma_f_235, sigma_a_235-sigma_f_235, sigma_a_238-sigma_f_238)
print(f'n = {thermal_fission_factor_start} at start')
#print(f'n = {thermal_fission_factor_end} at end of cycle')

resonance_passage_factor_start = calc_p(enrichment_percent, N_U, sigma_I_abs, xi_Sigma_s,3)
resonance_passage_factor_end = calc_p(enrichment_percent_end, N_U, sigma_I_abs, xi_Sigma_s,3)
print (f'p = {resonance_passage_factor_start} at start')
#print (f'p = {resonance_passage_factor_end} at end of cycle')

# Calculating macroscopic absorption cross section for moderator (H2O)
sigma_a_H = 0.332e-24    #Reactor Physics KSU
sigma_a_H2O = 2 * sigma_a_H     #Presuming absorption for Oxygen = 0 (Reactor Physics KSU)
Sigma_a_H2O = (rho_mod * NA) / Ma_mod * sigma_a_H2O

# Calculating macroscopic absorption cross section for fuel (UO2)
Sigma_a_fuel = sigma_a_235 * N_U235 + sigma_a_238 * N_U238  

# Calculating reactor poisoning from 135Xe. We assume all poisoning is due to Xenon 135.
sigma_a_Xe = 2.6e-18 
yield_Xe = 1.016e-3     #Independent fission yield at e = 0.0253 eV
yield_I = 6.437e-2      #Cumulative fission yield at e = 0.0253 eV. Because Te-135 has T1/2 = 19.6s we assume all Te-135 is I-135 at equilibrium. 
T_half_Xe =  33120
Lambda_Xe = math.log(2) / T_half_Xe

Sigma_f_fuel = sigma_f_235 * enrichment_percent * rho_UO2 * NA / Ma_UO2     #Assumption: All fission occurs from U-235.
neutron_flux = reactor_power_wh / (sigma_f_235 * energy_per_fission * N_U235_start)   
N_Xe = (yield_I + yield_Xe) * Sigma_f_fuel * neutron_flux / (Lambda_Xe + sigma_a_Xe * neutron_flux)

# Finally calculating Macroscopic absorption cross section and volume fraction for Xe. 
Sigma_a_Xe = N_Xe * sigma_a_Xe
V_rel_Xe = N_Xe / N_U   

# Next step is calculating the safety rods physical specifications. In Ringhals the safety rods are B4C.
rho_B4C = 2.52  #g/cm^3
Ma_B4C = 55.26  #g/mol
number_of_rods = 157    #From Ringhals specifications
blade_width_rod = 15   #cm Assumed values
thickness_rod = 0.5  #cm
N_B4C = rho_B4C*NA/Ma_B4C
sigma_a_bor = 3840e-24
Sigma_a_bor = N_B4C * sigma_a_bor


def calc_mod_and_rod_volumes(rod_insertion, scale_factor):
    length_rod = core_height_m * 100 * rod_insertion 

    #Relative and total volume of the control rods
    V_fuel = core_mass / rho_UO2    # Under the assumption that the full core is made up of fuel
    V_tot_rods = number_of_rods * (2 * length_rod * blade_width_rod  * thickness_rod)
    # scale_factor: Accounting for the self-shielding and decrease of neutron flux around the rod.

    # Calibrated by making sure k_eff = 1 at start. 
    V_rel_rods = V_tot_rods / V_fuel * scale_factor

    #Volume of the moderator. The difference between full volume of the core and the volumes of fuel and control rods is calculated. 
    V_core = math.pi * core_height_m * 100 * (core_diameter_m * 100 / 2)**2    # Actual volume of the core
    V_tot_mod = V_core - V_fuel - V_tot_rods    
    V_rel_mod = V_tot_mod / V_fuel
    #print(f'Volume mod: {V_tot_mod}')
    #print(f'Volume fuel: {V_fuel}')
    #print(f'Volume core: {V_core}')
    return V_rel_rods, V_rel_mod


# Calculation of the thermal utilization factor (f)
def calc_f(rod_insertion, scale_factor):
    V_rel_rods, V_rel_mod = calc_mod_and_rod_volumes(rod_insertion, scale_factor)
    
    f_fuel = Sigma_a_fuel 
    f_Xe = Sigma_a_Xe * V_rel_Xe
    f_mod = Sigma_a_H2O * V_rel_mod
    f_rods = Sigma_a_bor * V_rel_rods
    
    f = f_fuel / (f_fuel + f_mod + f_Xe + f_rods)
    #print("f_fuel:", f_fuel)
    #print("f_mod:", f_mod)
    #print("f_Xe:", f_Xe)
    #print("f_rods:", f_rods)
    
    return f


# ε (epsilon) - Fast fission factor (assumed reasonable level for UO2 in BWR)
epsilon = 1.06
print(f'ε = {epsilon}')

# P - Total non-leakage factor (estimate for large core, otherwise 1.0 if leakage is neglected)
P_leakage = 0.97
print (f'P = {P_leakage}')

# Function to calculate k_eff
def calc_k_eff(n, E, p, f, P):
    return n * E * p * f * P


# Making sure k = 1 at start and end by varying rod insertion. Also 
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
    
    # 4. Calculate volumes and (f) as the rod insertion changes
    V_rel_rods, V_rel_mod = calc_mod_and_rod_volumes(rod_insertion, scale_factor=0.15)
    f_val = sig_a_fuel / (sig_a_fuel + current_Sigma_a_H2O * V_rel_mod + Sigma_a_Xe * V_rel_Xe + Sigma_a_bor * V_rel_rods)
    
    return calc_k_eff(n_val, epsilon, p_val, f_val, P_leakage)

# Find rod insertion for k = 1
def target_critical(rod_insertion, current_enrichment):
    return evaluate_reactor_state(rod_insertion, current_enrichment) - 1.0

rod_critical_start = brentq(target_critical, 0.0, 0.99, args=(enrichment_percent,))
rod_critical_end = brentq(target_critical, 0.0, 0.99, args=(enrichment_percent_end,))

print("\n REACTOR STATE AT k=1 ")
print(f"Start (Enrichment {enrichment_percent*100:.2f}%): Control rods {rod_critical_start*100:.1f}% inserted")
print(f"End   (Enrichment {enrichment_percent_end*100:.2f}%): Control rods {rod_critical_end*100:.1f}% inserted")


# Ploting void coefficients for the start and end of the cycle to prove that increasing the void would make the reactor undercritical. 

void_fractions = np.linspace(0.0, 0.8, 50)  # Test void from 0% to 80%
k_eff_start_list = []
k_eff_end_list = []

for v in void_fractions:
    # Call reactor function with v as void so that it varies. 
    k_start = evaluate_reactor_state(rod_critical_start, enrichment_percent, current_void=v)
    k_end = evaluate_reactor_state(rod_critical_end, enrichment_percent_end, current_void=v)
    
    k_eff_start_list.append(k_start)
    k_eff_end_list.append(k_end)

# Plotting the results
plt.figure(figsize=(10, 6))
plt.plot(void_fractions * 100, k_eff_start_list, label="Start of Cycle (BOC)", color="blue", linewidth=2)
plt.plot(void_fractions * 100, k_eff_end_list, label="End of Cycle (EOC)", color="red", linewidth=2, linestyle='--')

# Mark void point (40%), which is the void we have assumed in the assigment. 
plt.axvline(40, color='gray', linestyle=':', label='Operational Void (40%)')
plt.axhline(1.0, color='black', linestyle='-', linewidth=1, label='Criticality (k = 1)')

plt.title("Void Coefficient: $k_{eff}$ vs Void Fraction", fontsize=14)
plt.xlabel("Void Fraction (%)", fontsize=12)
plt.ylabel("Effective Multiplication Factor ($k_{eff}$)", fontsize=12)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.show()