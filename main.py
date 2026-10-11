"""The task is to make a preliminary design of a reactor so that it is critical (k=1) when it has reached equilibrium at full power at start-up and can be operated until the next fuel change. In the design work, you must choose the following parameters: reactor power, reactor size, operating time between fuel changes, reactor type, fuel, enrichment, moderator, moderator temperature, void, fuel dimension and moderator/fuel volume ratio. Explain how you use control rods, burnable absorbers and/or boron. The moderator temperature coefficient / void coefficient should be negative both at start-up and at the end of the operating period. Some of your own assumptions will need to be made, especially regarding neutron absorption in reactor poisons and construction materials."""
import math as math
import matplotlib.pyplot as plt
import numpy as np

NA = 6.022 * 10 ** 23  # Avogadros tal, st/mol

#Technical specifications for Ringhals 1
reactor_power_mwth = 2540.0         
core_assemblies = 648
assemblies_mass = 193e3                
core_height_m = 3.7                 
core_diameter_m = 5.95             
core_mass = core_assemblies * assemblies_mass #kg
op_time_between_refueling_mo = 12   
reactor_type = "BWR"  


#fuel specifications
fuel_material = "UO2"  
Ma_UO2 = 270  # g/mol UO2
rho_UO2 = 10.97  # g/cm^3             
enrichment_percent = 3.8/100
N_U = (rho_UO2 * NA) / Ma_UO2
N_U235 = enrichment_percent * N_U
N_U238 = (1 - enrichment_percent) * N_U 

#Fissile and fertile cores at the start of the cycle
energy_per_fission = 3.2e11     # Joules / Fission
energy_per_cycle = reactor_power_mwth * 60 * 60 * 24 * 365  # Assumption: 1 cycle = 1 year
fissions_per_cycle = energy_per_cycle / energy_per_fission 
N_U235_start = enrichment_percent * core_mass * NA / Ma_UO2  #Number of fissile cores at the start of the cycle. 
                                                             #Assumption: Only fissile cores are U-235
N_U238_start = (1 - enrichment_percent) * core_mass * NA / Ma_UO2

#moderator specifications
moderator = "H2O"                  
moderator_temp_c = (286.0 + 273)    
void_fraction = 0.40                
rho_steam = 36.53e-3  # g/cm^3 at T = 286 C and P = 7 MPa from
rho_water = 739.7e-3  # g/cm^3
rho_mod = void_fraction * rho_steam + (1 - void_fraction) * rho_water
Ma_mod = 18  # g/cm^3 för H2O
fuel_pin_diameter_mm = 9.8          
mod_fuel_volume_ratio = 3        #(Vm/Vf) assumed given lectures in course
sigma_I_abs = 20.0 * 10**-24    # barn converted to cm^2
xi_Sigma_s = 0.8                # cm^-1 (assumed value for H2O)

# cross sections from nndc.bnl.gov
sigma_f_235 = 596 * 10 ** -24  # barn
sigma_a_235 = 97.6 * 10 ** -24 + sigma_f_235  #(n,gamma) + (n, total fission)  # barn
sigma_f_238 = 1.76 * 10 ** -5 * 10 ** -24
sigma_a_238 = 2.78 * 10 ** -24 + sigma_f_238


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

thermal_fission_factor = calc_n(2.43, enrichment_percent, sigma_f_235, sigma_a_235-sigma_f_235, sigma_a_238-sigma_f_238)
print(f'n = {thermal_fission_factor}')


resonance_passage_factor = calc_p(enrichment_percent, N_U, sigma_I_abs, xi_Sigma_s,3)
print (f'p = {resonance_passage_factor}')

# Calculating macroscopic absorption cross section for moderator (H2O)
sigma_a_H = 0.332e-24    #Reaktorfysik KSU
sigma_a_H2O = 2 * sigma_a_H     #Presuming absorption for Oxygen = 0 (Reaktorfysik KSU)
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
neutron_flux = reactor_power_mwth / (sigma_f_235 * energy_per_fission * N_U235_start)   
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
    #print(f'Volym mod: {V_tot_mod}')
    #print(f'Volym fuel: {V_fuel}')
    #print(f'Volym core: {V_core}')
    return V_rel_rods, V_rel_mod


# Beräkning av termiska utnyttjandefaktorn (f)
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

rod_insertion_start = 0.0
f_val = calc_f(rod_insertion_start, scale_factor = 0.15)
print(f"f vid {rod_insertion_start*100}% inskjutning: {f_val:.4f}")

# Beräkna de fem faktorerna

# ε (epsilon) - Snabbfissionsfaktorn (antagen rimlig nivå för UO2 i BWR)
epsilon = 1.02
print(f'ε = {epsilon}')

# P - Total icke-läckagefaktor (uppskattning för stor härd, annars 1.0 om läckage försummas)
P_leakage = 0.97
print (f'P = {P_leakage}')

# Funktion för att beräkna k_eff
def calc_k_eff(n, E, p, f, P):
    return n * E * p * f * P

# Kör beräkning för 0% inskjutning (helt utdragna stavar) 
print("\n--- Körning vid 0% inskjutning ---")
f_val_0 = calc_f(rod_insertion=0.0, scale_factor=0.1254)
k_eff_0 = calc_k_eff(thermal_fission_factor, epsilon, resonance_passage_factor, f_val_0, P_leakage)
print(f"f vid 0%: {f_val_0:.4f}")
print(f"k_eff vid 0% inskjutning: {k_eff_0:.4f}")

# --- 4. Kör beräkning för 70% inskjutning (driftstart) ---
print("\n--- Körning vid 70% inskjutning ---")
f_val_70 = calc_f(rod_insertion=0.70, scale_factor=0.1254)
k_eff_70 = calc_k_eff(thermal_fission_factor, epsilon, resonance_passage_factor, f_val_70, P_leakage)
print(f"f vid 70%: {f_val_70:.4f}")
print(f"k_eff vid 70% inskjutning: {k_eff_70:.4f}")