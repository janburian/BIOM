import math

def Nernst_equation(z, X_ext, X_int, temperature):
    # Constants
    R = 8.31447  # J/(mol·K)
    F = 96485.3321  # C/mol

    T_kelvin = temperature + 273.15

    E_equilibrium = (R * T_kelvin) / (z * F) * math.log(X_ext / X_int)  # [mV]

    return E_equilibrium * 1_000  # [V]


temperature = 37  # Temperature in Celsius

# Sodium
z_Na = 1
Na_ext = 142
Na_int = 15

# Potassium
z_K = 1
K_ext = 4
K_int = 150

# Leakage (Cl)
z_L = -1
L_ext = 120
L_int = 5


E_Na = Nernst_equation(z_Na, Na_ext, Na_int, temperature)
E_K = Nernst_equation(z_K, K_ext, K_int, temperature)
E_L = Nernst_equation(z_L, L_ext, L_int, temperature)

print(f"Equilibrium potential for sodium: E_Na = {E_Na} mV")
print(f"Equilibrium potential for potassium: E_K = {E_K} mV")
print(f"Equilibrium potential for leakage: E_L = {E_L} mV")