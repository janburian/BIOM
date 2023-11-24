# Import modules
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# Hodgkin-Huxley model parameters
C_m = 1.0  # membrane capacitance (uF/cm^2)

g_Na = 120.0  # sodium conductance (mS/cm^2)
g_K = 36.0  # potassium conductance (mS/cm^2)
g_L = 0.3  # leak conductance (mS/cm^2)

E_Na = 60.0  # sodium reversal potential (mV)
E_K = -77.0  # potassium reversal potential (mV)
E_L = -54.0  # leak reversal potential (mV)


# Function representing the Hodgkin-Huxley model
def hodgkin_huxley(y, t):
    V_m, m, h, n = y

    # Membrane currents
    I_Na = g_Na * m**3 * h * (V_m - E_Na)
    I_K = g_K * n**4 * (V_m - E_K)
    I_L = g_L * (V_m - E_L)

    # Total membrane current
    I_ion = I_Na + I_K + I_L

    # Differential equations
    dV_dt = (1 / C_m) * (I_inj(t) - I_ion)

    dn_dt = alpha_n(V_m) * (1 - n) - beta_n(V_m) * n
    dm_dt = alpha_m(V_m) * (1 - m) - beta_m(V_m) * m
    dh_dt = alpha_h(V_m) * (1 - h) - beta_h(V_m) * h

    return [dV_dt, dn_dt, dm_dt, dh_dt]


# External current injection function
def I_inj(t):
    # You can modify this function to simulate different current injections
    return 20.0 if 5 < t < 15 else 0.0


# Alpha and beta functions for m, h, and n gates
def alpha_m(V):
    return 0.1 * (V + 40.0) / (1.0 - np.exp(-(V + 40.0) / 10.0))


def beta_m(V):
    return 4.0 * np.exp(-(V + 65.0) / 18.0)


def alpha_h(V):
    return 0.07 * np.exp(-(V + 65.0) / 20.0)


def beta_h(V):
    return 1.0 / (1.0 + np.exp(-(V + 35.0) / 10.0))


def alpha_n(V):
    return 0.01 * (V + 55.0) / (1.0 - np.exp(-(V + 55.0) / 10.0))


def beta_n(V):
    return 0.125 * np.exp(-(V + 65.0) / 80.0)


if __name__ == "__main__":
    # Determine the initial conditions
    resting_potential = -65.0  # mV
    initial_m = 0.05
    initial_h = 0.6
    initial_n = 0.3

    initial_conditions = [resting_potential, initial_m, initial_h, initial_n]

    # Time vector (Time of simulation)
    time = np.arange(0, 50, 0.01)

    # Solve the differential equations using odeint
    solution = odeint(hodgkin_huxley, initial_conditions, time)

    # Plot the results
    time_inj = np.arange(0, 50, 0.01)
    injection_values = [I_inj(t) for t in time_inj]

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.plot(time_inj, injection_values, label='I_inj')
    plt.title('External Current Injection')
    plt.xlabel('Time (ms)')
    plt.ylabel('Current (uA/cm^2)')
    plt.legend()


    plt.figure(figsize=(10, 6))
    plt.plot(time, solution[:, 0], label='V')
    plt.title('Hodgkin-Huxley Model')
    plt.xlabel('Time (ms)')
    plt.ylabel('Membrane Potential (mV)')
    plt.legend()
    plt.show()