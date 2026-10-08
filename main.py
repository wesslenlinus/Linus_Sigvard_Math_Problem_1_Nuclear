"""The task is to make a preliminary design of a reactor so that it is critical (k=1) when it has reached equilibrium at full power at start-up and can be operated until the next fuel change. In the design work, you must choose the following parameters: reactor power, reactor size, operating time between fuel changes, reactor type, fuel, enrichment, moderator, moderator temperature, void, fuel dimension and moderator/fuel volume ratio. Explain how you use control rods, burnable absorbers and/or boron. The moderator temperature coefficient / void coefficient should be negative both at start-up and at the end of the operating period. Some of your own assumptions will need to be made, especially regarding neutron absorption in reactor poisons and construction materials."""

#Technical specifications for Ringhals 1 and assumed values
reactor_power_mwth = 2540.0         # Termisk reaktoreffekt i MW
core_assemblies = 648               # Reaktorstorlek: Antal bränslepatroner
core_height_m = 3.7                 # Reaktorstorlek: Aktiv härdhöjd i meter
core_diameter_m = 4.2               # Reaktorstorlek: Härddiameter i meter
op_time_between_refueling_mo = 12   # Drifttid mellan bränslebyten i månader
reactor_type = "BWR"                # Reaktortyp
fuel_material = "UO2"               # Bränsle
enrichment_percent = 3.8            # Anrikning U-235 i procent (snitt)
moderator = "H2O"                   # Moderator
moderator_temp_c = 286.0            # Moderatortemperatur i Celsius
void_fraction = 0.40                # Genomsnittlig void (ånghalt)
fuel_pin_diameter_mm = 9.8          # Bränsledimension: stavdiameter i mm
mod_fuel_volume_ratio = 2.25        # Volymförhållande moderator/bränsle (Vm/Vf)



def k(n, E, p, f, P):
    k = n*E*p*f*P
    return (k)

# Beräkna mikroskopiskt absorptionstvärsnitt för moderator (H2O)
sigma_a_H = 0.332e-24    #Hämtat ur Reaktorfysik KSU
sigma_a_H20 = 2 * sigma_a_H     #Antar att syre har mikroskopiskt absorptionstvärsnitt = 0 enligt Reaktorfysik KSU

sigma_a_br =    # Mikroskopiskt absorptionstvärsnitt för bränsle

Sigma_a_br =    # Makroskopiskt absorptionstvärsnitt för bränsle

def f(S_a_br, V_br, S_a_mod, V_a_mod, S_a_styrstav, V_a_styrstav, S_a_Xe, V_a_Xe):
    