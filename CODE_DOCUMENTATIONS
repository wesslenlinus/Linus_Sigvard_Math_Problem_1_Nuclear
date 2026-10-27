## Assumptions (same as in LaTeX report)
The following assumptions have been made:

### Reactor Parameters (Ringhals 1)
* **Reactor Type:** BWR
* **Thermal Power:** 2540 MW
* **Core Dimensions (H x D):** 3.7 m x 5.95 m
* **Number of Fuel Assemblies:** 648
* **Operating Cycle:** 12 months

### Fuel and Control Rod Parameters
* **Fuel Material:** UO₂
* **Enrichment (²³⁵U):** 3.0 %
* **Fuel Density (ρ_UO₂):** 10.97 g/cm³
* **Control Rod Material:** B₄C
* **Control Rod Density (ρ_B₄C):** 2.52 g/cm³
* **Number of Control Rods:** 157
* **Calibration Factor:** 0.15

### Moderator Parameters
* **Moderator Material:** H₂O (light water)
* **Temperature:** 286 °C
* **Void Fraction:** 40 %
* **Volume Ratio (V_m/V_u):** 3

### Physical Cross Sections and Constants
* **Energy per Fission:** 3.2 × 10⁻¹¹ J
* **Effective Resonance Integral (σ_{I,abs}):** 20 barn
* **Slowing-down Power (ξΣ_s):** 0.8 cm⁻¹
* **Fission Cross Section (σ_{f,235}):** 596 barn
* **Absorption H₂O (σ_{a,H2O}):** 0.664 barn
* **Cumulative Iodine Yield (γ_I):** 6.437 %
* **Independent Xenon Yield (γ_{Xe}):** 0.1016 %

## Modeling Assumptions for the Five-Factor Formula (same as in LaTeX report)

### General Reactor & Fuel Usage
* **Sole Fissile Isotope:** Uranium-235 is assumed to be the only fissile isotope throughout the entire operating cycle. The model explicitly neglects the breeding and subsequent fissioning of Plutonium-239.
* **Power Generation:** All thermal power is assumed to be generated exclusively by the fission of U-235, neglecting minor energy contributions from U-238 fast fissions.
* **Neutrons per Fission ($\bar{\nu}$):** Treated as a constant value of 2.43 (representing pure U-235).

### Multiplication Factors
* **Fast Fission Factor ($\epsilon$):** Assumed to be a constant value of 1.06.
* **Non-Leakage Probability ($P$):** Assumed to be a constant value of 0.97.

### Resonance Escape Probability ($p$)
* **Effective Resonance Integral ($\sigma_{I,abs}$):** Assumed to be a constant 20 barn.
* **Slowing-Down Power ($\xi\Sigma_s$):** Assumed to be a constant 0.8 cm⁻¹. This static value implicitly ignores dynamic local variations in moderator density caused by the 40% void fraction.

### Thermal Utilization Factor ($f$)
* **Active Absorbers:** Thermal neutron capture is assumed to occur *only* in the fuel (U-235, U-238), the moderator (H₂O), Xenon-135 (treated as the sole reactor poison), and the B₄C control rods.
* **Ignored Materials:** Absorption from Zircaloy cladding, internal structural components, and other fission poisons (such as Samarium-149) is strictly neglected.
* **Control Rod Self-Shielding:** Because B₄C is an exceptionally strong absorber, an empirical scaling factor of 0.15 is applied to the physical control rod volume to account for spatial self-shielding (flux depression) at the blade surfaces.
* **Moderator Volume:** The available water moderator is calculated dynamically as the remaining empty space in the core after subtracting the fixed fuel volume and the currently inserted control rod volume.
