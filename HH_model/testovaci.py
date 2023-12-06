import numpy as np
import matplotlib.pyplot as plt

# Constants
C_m = 1.0
g_Na = 120.0
g_K = 36.0
g_L = 0.3
E_Na = 50.0
E_K = -77.0
E_L = -54.387

# Simulation parameters
dt = 0.01
t_max = 50.0
time_points = np.arange(0, t_max, dt)

# Alpha and beta functions for m, h, and n gates
def alpha_n(V):
    alpha_n = (0.01 * (V + 10)) / (np.exp((V + 10) / 10) - 1)
    return alpha_n


def beta_n(V):
    beta_n = 0.125 * np.exp(V / 80)
    return beta_n


def alpha_m(V):
    alpha_m = 0.1 * (V + 25) / (np.exp((V + 25) / 10) - 1)
    return alpha_m


def beta_m(V):
    beta_m = 4 * np.exp(V / 18)
    return beta_m


def alpha_h(V):
    alpha_h = 0.07 * np.exp(V / 20)
    return alpha_h


def beta_h(V):
    beta_h = 1 / (np.exp((V + 30) / 10) + 1)
    return beta_h

# Initial conditions
V = -65.0
m = alpha_m(V) / (alpha_m(V) + beta_m(V))
h = alpha_h(V) / (alpha_h(V) + beta_h(V))
n = alpha_n(V) / (alpha_n(V) + beta_n(V))

# Simulation loop
V_trace = []

for t in time_points:
    # Update gating variables using alpha and beta functions
    m += (alpha_m(V) * (1 - m) - beta_m(V) * m) * dt
    h += (alpha_h(V) * (1 - h) - beta_h(V) * h) * dt
    n += (alpha_n(V) * (1 - n) - beta_n(V) * n) * dt

    # Compute membrane current based on gating variables
    I_Na = g_Na * m**3 * h * (V - E_Na)
    I_K = g_K * n**4 * (V - E_K)
    I_L = g_L * (V - E_L)

    # Total membrane current
    I_total = I_Na + I_K + I_L

    # Update membrane potential using the current
    V += (I_total / C_m) * dt
    V_trace.append(V)

# Plot the results
plt.figure(figsize=(10, 6))
plt.plot(time_points, V_trace, label='V')
plt.title('Hodgkin-Huxley Model Simulation')
plt.xlabel('Time (ms)')
plt.ylabel('Membrane Potential (mV)')
plt.legend()
plt.show()