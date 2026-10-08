"""The task is to make a preliminary design of a reactor so that it is critical (k=1) when it has reached equilibrium at full power at start-up and can be operated until the next fuel change. In the design work, you must choose the following parameters: reactor power, reactor size, operating time between fuel changes, reactor type, fuel, enrichment, moderator, moderator temperature, void, fuel dimension and moderator/fuel volume ratio. Explain how you use control rods, burnable absorbers and/or boron. The moderator temperature coefficient / void coefficient should be negative both at start-up and at the end of the operating period. Some of your own assumptions will need to be made, especially regarding neutron absorption in reactor poisons and construction materials."""

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
enrichment_percent = 3.8

#moderator specifications
moderator = "H2O"                  
moderator_temp_c = (286.0 + 273)    
void_fraction = 0.40                
fuel_pin_diameter_mm = 9.8          
mod_fuel_volume_ratio = 3        #(Vm/Vf) assumed given



def k(n, E, p, f, P):
    k = n*E*p*f*P
    return (k)

# cross sections from nndc.bnl.gov
sigma_f_235 = 596 * 10 ** -24  # barn
sigma_a_235 = 97.6 * 10 ** -24 + sigma_f_235  #(n,gamma) + (n, total fission)  # barn
sigma_f_238 = 1.76 * 10 ** -5 * 10 ** -24
sigma_a_238 = 2.78 * 10 ** -24 + sigma_f_238

#calculation of thermal fission factor n
def n(newbar, N235, N238, sigmaF235, sigmaC235, sigmaC238):
    return(n)