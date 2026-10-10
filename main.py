"""The task is to make a preliminary design of a reactor so that it is critical (k=1) when it has reached equilibrium at full power at start-up and can be operated until the next fuel change. In the design work, you must choose the following parameters: reactor power, reactor size, operating time between fuel changes, reactor type, fuel, enrichment, moderator, moderator temperature, void, fuel dimension and moderator/fuel volume ratio. Explain how you use control rods, burnable absorbers and/or boron. The moderator temperature coefficient / void coefficient should be negative both at start-up and at the end of the operating period. Some of your own assumptions will need to be made, especially regarding neutron absorption in reactor poisons and construction materials."""
import math as math

NA = 6.022 * 10 ** 23  # Avogadros tal, st/mol

#Technical specifications for Ringhals 1
reactor_power_mwth = 2540.0         
core_assemblies = 648               
core_height_m = 3.7                 
core_diameter_m = 4.2               
op_time_between_refueling_mo = 12   
reactor_type = "BWR"  

#fuel specifications
fuel_material = "UO2"  
Ma_UO2 = 270  # g/mol UO2
rho_UO2 = 10.97  # g/cm^3             
enrichment_percent = 3.8/100

#moderator specifications
moderator = "H2O"                  
moderator_temp_c = (286.0 + 273)    
void_fraction = 0.40                
rho_steam = 36.53 / 1e3  # g/cm^3 at T = 286 C and P = 7 MPa from
rho_water = 739.7 / 1e3  # g/cm^3
rho_mod = void_fraction * rho_steam + (1 - void_fraction) * rho_water
Ma_mod = 18  # g/cm^3 för H2O
fuel_pin_diameter_mm = 9.8          
mod_fuel_volume_ratio = 3        #(Vm/Vf) assumed given lectures in course
N_U = (rho_UO2 * NA) / Ma_UO2
N_U235 = enrichment_percent * N_U
N_U238 = (1 - enrichment_percent) * N_U 
sigma_I_abs = 20.0 * 10**-24    # barn converted to cm^2
xi_Sigma_s = 0.8                # cm^-1 (assumed value for H2O)


def k(n, E, p, f, P):
    k = n*E*p*f*P
    return (k)

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
print(thermal_fission_factor)


resonance_passage_factor = calc_p(enrichment_percent, N_U, sigma_I_abs, xi_Sigma_s,3)
print (resonance_passage_factor)

# Calculating macroscopic absorption cross section for moderator (H2O)
sigma_a_H = 0.332e-24    #Reaktorfysik KSU
sigma_a_H2O = 2 * sigma_a_H     #Presuming absorption for Oxygen = 0 (Reaktorfysik KSU)
Sigma_a_H2O = (rho_mod*NA)/Ma_mod * sigma_a_H2O

# Calculating macroscopic absorption cross section for fuel (UO2)
Sigma_a_br = sigma_a_235*N_U235 + sigma_a_238*N_U238  



def f(S_a_br, V_br, S_a_mod, V_a_mod, S_a_styrstav, V_a_styrstav, S_a_Xe, V_a_Xe):
    