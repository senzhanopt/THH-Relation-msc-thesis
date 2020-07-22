# This code aims to find the T-H-H relationship of a typical alkaline electrolyser
# with 30% mass KOH at a pressure of 30 bar
from numpy import exp, log
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 2                   # electrons transferred
F = 96485               # Faraday constant
p = 30 * 0.986923267    # pressure [atm]
R = 8.314               # gas constant
M = 8                   # % molarity [mol/l]
ilim = 30               # limiting current density [A/cm^2]

def voltage_thermalNeutral(t): # Celsius degree
    T = 273.15 + t # [K]
    voltage_HHV = 1.4756 + 2.252E-4 * t + 1.52E-8 * t**2  # (V)
    Y = 42960 + 40.762 * t - 0.06682 * t**2  # (J/mol)
    m = 7.64 # molality for 30% KOH solution (mol/kg)
    pw_water = exp(37.04 - 6276/T - 3.416 * log(T)) #vapour pressure of water
    pw_KOH = exp(0.01621 - 0.1380 * m + 0.1933 * m**0.5 + 1.024 * log(pw_water)) #vapour pressure of KOH solution
    Phi = 1.5 * pw_KOH / (p - pw_KOH)
    return voltage_HHV + Phi * Y / (n * F)

def voltage_cell(t, i): # [A/cm^2]
    T = 273.15 + t # [K]
    m = 7.64 # molality for 30% KOH solution (mol/kg)
    pw_water = exp(37.04 - 6276/T - 3.416 * log(T)) #vapour pressure of water
    pw_KOH = exp(0.01621 - 0.1380 * m + 0.1933 * m**0.5 + 1.024 * log(pw_water)) #vapour pressure of KOH solution
    voltage_rev = 1.5184 - 1.5421E-3 * T + 9.523E-5 * T * log(T) + 9.84E-8 * T**2
    voltage_rev += R * T / (n * F) * log((p-pw_KOH)**1.5 * pw_water / pw_KOH)
    alpha_anode = 0.0675 + 0.00095*T # anode (Ni) charge-transfer coefficient
    alpha_cathode = 0.1175 + 0.00095*T # cathode (Ni) charge-transfer coefficient
    i0 = 3.15e-4 # exchange current density [A/cm^2]
    Vact_anode = R*T/(n*F*alpha_anode)*log(i/i0)
    Vact_cathode = R*T/(n*F*alpha_cathode)*log(i/i0)
    delta = 0.66 #  electrolyte thickness [cm]
    sigma = -2.041 * M - 0.0028*M**2 + 0.001043*M**3 + 0.005332*M*T + 207.2*M/T - 0.0000003*M**2*T**2 # conductivity of KOH solution [S/cm]
    # Calculate the effect of bubble phenomena
    eps = 0.0153*(i/ilim)**0.3 # void fraction of the electrolyte, porosity
    sigma_eps = (1-eps)**1.5*sigma # electrical conductivity in the presence of bubbles [S/cm]
    r = delta/sigma_eps # area specific resistance [cm^2/S]
    Vohm = i * r
    return voltage_rev + Vact_anode + Vact_cathode + Vohm

t_set = []
P_H2_set = []
P_heat_set = []
cellArea = 2500 # [cm^2]
for t in np.linspace(60, 80, 21):
    for i in np.linspace(0.2, 0.4, 21):
        t_set.append(t)
        P_H2_set.append(voltage_thermalNeutral(t)*i*cellArea / 1000)
        P_heat_set.append((voltage_cell(t,i)-voltage_thermalNeutral(t))*i*cellArea / 1000)

fig = plt.figure(figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(t_set, P_H2_set, P_heat_set, c='b', marker='*')
plt.grid()
ax.text(t_set[0], P_H2_set[0], P_heat_set[0], 'A', fontsize = 20)
ax.text(t_set[20], P_H2_set[20], P_heat_set[20], 'B', fontsize = 20)
ax.text(t_set[440], P_H2_set[440], P_heat_set[440], 'C', fontsize = 20)
ax.text(t_set[420], P_H2_set[420], P_heat_set[420], 'D', fontsize = 20)
ax.set_xlabel('t[$^\circ$C]', fontsize = 18)
ax.set_ylabel('$P_{H2}$[kW]', fontsize = 18)
ax.set_zlabel('$P_{Heat}$[kW]', fontsize = 18)
plt.savefig("thh.eps")
plt.show()

fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
for t in [60,65,70,75,80]:
    ax2.plot(P_H2_set[(t-60)*21:(t-59)*21], P_heat_set[(t-60)*21:(t-59)*21], label = "{}$^\circ$C".format(t))
ax2.set_xlabel('$P_{H2}$[kW]', fontsize = 12)
ax2.set_ylabel('$P_{Heat}$[kW]', fontsize = 12)
ax2.legend()
plt.savefig("hh.eps")
plt.show()

fig3 = plt.figure(figsize=(10,7))
ax3 = fig3.add_subplot(111, projection='3d')
ax3.plot([t_set[0], t_set[20],t_set[440],t_set[420],t_set[0]],[P_H2_set[0], P_H2_set[20],P_H2_set[440],P_H2_set[420],P_H2_set[0]],[P_heat_set[0], P_heat_set[20],P_heat_set[440],P_heat_set[420],P_heat_set[0]], c='b')
ax3.plot([t_set[0], t_set[220]],[P_H2_set[0], P_H2_set[220]],[P_heat_set[0], P_heat_set[220]], 'b')
ax3.plot([t_set[20], t_set[220]],[P_H2_set[20], P_H2_set[220]],[P_heat_set[20], P_heat_set[220]], 'b')
ax3.plot([t_set[440], t_set[220]],[P_H2_set[440], P_H2_set[220]],[P_heat_set[440], P_heat_set[220]], 'b')
ax3.plot([t_set[420], t_set[220]],[P_H2_set[420], P_H2_set[220]],[P_heat_set[420], P_heat_set[220]], 'b')
ax3.text(t_set[0], P_H2_set[0], P_heat_set[0], 'A', fontsize = 20)
ax3.text(t_set[20], P_H2_set[20], P_heat_set[20], 'B', fontsize = 20)
ax3.text(t_set[440], P_H2_set[440], P_heat_set[440], 'C', fontsize = 20)
ax3.text(t_set[420], P_H2_set[420], P_heat_set[420], 'D', fontsize = 20)
ax3.text(t_set[220], P_H2_set[220], P_heat_set[220], 'E', fontsize = 20)
ax3.set_xlabel('t[$^\circ$C]', fontsize = 18)
ax3.set_ylabel('$P_{H2}$[kW]', fontsize = 18)
ax3.set_zlabel('$P_{Heat}$[kW]', fontsize = 18)
plt.savefig("thh2.eps")
plt.show()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
for t in [60,65,70,75,80]:
    ax4.plot(np.linspace(0.2,0.4,21), [(P_heat_set[ind]/P_H2_set[ind]) for ind in range((t-60)*21, (t-59)*21)], label = "{}$^\circ$C".format(t))
ax4.set_xlabel('current density[$A/cm^2$]', fontsize = 12)
ax4.set_ylabel('$P_{Heat}/P_{H2}$', fontsize = 12)
ax4.legend()
plt.savefig("ratio.eps")
plt.show()
