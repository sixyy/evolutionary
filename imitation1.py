import math
# from geomdl.visualization import VisMPL as vis
from geomdl import fitting
# import pandas as pd    #Processing Data from Spreadsheet
import numpy as np
# from scipy.integrate import odeint #Ordinary differential equation solver
import matplotlib.pyplot as plt    #For graph plotting
# from matplotlib.patches import FancyArrowPatch
from sympy import *
from scipy.integrate import odeint
from scipy.interpolate import make_interp_spline

x, n = symbols('x n')

# fig2-1
R0 = 6.0
S0 = 6.0
T0 = 4.0
P0 = 1.0
R1 = 2.0
S1 = 4.0
T1 = 3.0
P1 = 6.0
theta = 0.1
episi = 1.0
beta = 1.0
y0 = [0.9, 0.5]
y1 = [0.01, 0.01]

# fig3-1
# R0 = 2.0
# S0 = 9.0
# T0 = 1.0
# P0 = 7.0
# R1 = 6.0
# S1 = 1.0
# T1 = 9.0
# P1 = 9.0
# theta = 1.0
# episi = 1.0
# beta = 1.0
# y0 = [0.5, 0.99]
# y1 = [0.5, 0.18]

t = np.linspace(0, 1000, 10000)

def pend(y, t, R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta):
    x, n = y
    dydt = [x*(1 - x)*(1/(1 + math.exp(-beta* ((((R1 - R0 - S1 + S0)*n + (R0 - S0))*x + (S1 - S0)*n + S0) - (((T1 - T0 - P1 + P0)*n + (T0 - P0))* x + (P1 - P0)*n + P0)))) - 1/(1 + math.exp(beta *((((R1 - R0 - S1 + S0)*n + (R0 - S0))*x + (S1 - S0)*n + S0) - (((T1 - T0 - P1 + P0)*n + (T0 - P0))*x + (P1 - P0)*n + P0))))), episi*n*(1 - n)*(-1 + (1 + theta)*x)]
    return dydt


sol = odeint(pend, y0, t, args=(R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta))
sol1 = odeint(pend, y1, t, args=(R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta))

font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 35}

fig1 = plt.figure(figsize=(10, 10))
bwith = 2 #边框宽度设置为2
ax = plt.gca() #获取边框
ax.spines['bottom'].set_linewidth(bwith)
ax.spines['left'].set_linewidth(bwith)
ax.spines['top'].set_linewidth(bwith)
ax.spines['right'].set_linewidth(bwith)


plt.plot(sol[:, 0], sol[:, 1], 'lightgrey')
plt.plot(sol1[:, 0], sol1[:, 1], 'lightgrey')
plt.scatter(y0[0], y0[1], c='dodgerblue', s=200, edgecolors='dodgerblue', linewidths=2, zorder = 10)
plt.scatter(y1[0], y1[1], c='red', s=200, edgecolors='red', linewidths=2, zorder = 10)
plt.scatter(1/(1 + theta), ((T0 - R0) + theta * (P0 - S0)) / ((R1 - T1 + T0 - R0) + theta * (S1 - P1 + P0 - S0)), c='black', s=200, edgecolors='black', linewidths=2, zorder = 10)
# print(1/(1 + theta), ((T0 - R0) + theta * (P0 - S0)) / ((R1 - T1 + T0 - R0) + theta * (S1 - P1 + P0 - S0)))

# delta = 0.01
# xrange = np.arange(0, 1, delta)
# yrange = np.arange(0, 1, delta)
# x, n = np.meshgrid(xrange,yrange)
# plt.contour(x, n, -R0*(x - (a + 1)/(a + 2))*(c*(a + 1)*(VH*a + 2*VH - VL*a - 2*VL) - sqrt(-(VH - VL)*(a + 1)*(a + 2)*(4*c**2*R0**2 - 8*c**2*R0*a - 16*c**2*R0 - c**2*VH*a**2 - 3*c**2*VH*a - 2*c**2*VH + c**2*VL*a**2 + 3*c**2*VL*a + 2*c**2*VL + 4*c**2*a**2 + 16*c**2*a + 16*c**2 - 4*c*R0**2*VH - 4*c*R0**2*VL + 4*c*R0*VH*a + 8*c*R0*VH + 4*c*R0*VL*a + 8*c*R0*VL + 4*R0**2*VH*VL)))*(VH**2 - 2*VH*VL + VL**2)/(2*(VH*a + 2*VH - VL*a - 2*VL)*(c**2*R0**2*a + 2*c**2*R0**2 - 2*c**2*R0*a**2 - 8*c**2*R0*a - 8*c**2*R0 + c**2*a**3 + 6*c**2*a**2 + 12*c**2*a + 8*c**2 - c*R0**2*VH*a - 2*c*R0**2*VH - c*R0**2*VL*a - 2*c*R0**2*VL + c*R0*VH*a**2 + 4*c*R0*VH*a + 4*c*R0*VH + c*R0*VL*a**2 + 4*c*R0*VL*a + 4*c*R0*VL + R0**2*VH*VL*a + 2*R0**2*VH*VL)) + n - (VH - VL)*(n - (-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)))*(a + 1)*(0.5*x - 0.5*(a + 1)/(a + 2))/(a + 2)**2 + (x - (a + 1)/(a + 2))*(-c*(a + 1)*(0.5*x - 0.5*(a + 1)/(a + 2))/R0 + (0.5*n - 0.5*(-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)))*(a + 2)*(c*(-R0 + a + 2) + R0*VH)*(c*(-R0 + a + 2) + R0*VL)/(R0**2*(VH - VL)**2)) - (-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)),0, colors='red', linewidths=2, linestyles='dashdot')
# # plt.contour(x, n, -R0*(x - (a + 1)/(a + 2))*(c*(a + 1)*(VH*a + 2*VH - VL*a - 2*VL) + sqrt(-(VH - VL)*(a + 1)*(a + 2)*(4*c**2*R0**2 - 8*c**2*R0*a - 16*c**2*R0 - c**2*VH*a**2 - 3*c**2*VH*a - 2*c**2*VH + c**2*VL*a**2 + 3*c**2*VL*a + 2*c**2*VL + 4*c**2*a**2 + 16*c**2*a + 16*c**2 - 4*c*R0**2*VH - 4*c*R0**2*VL + 4*c*R0*VH*a + 8*c*R0*VH + 4*c*R0*VL*a + 8*c*R0*VL + 4*R0**2*VH*VL)))*(VH**2 - 2*VH*VL + VL**2)/(2*(VH*a + 2*VH - VL*a - 2*VL)*(c**2*R0**2*a + 2*c**2*R0**2 - 2*c**2*R0*a**2 - 8*c**2*R0*a - 8*c**2*R0 + c**2*a**3 + 6*c**2*a**2 + 12*c**2*a + 8*c**2 - c*R0**2*VH*a - 2*c*R0**2*VH - c*R0**2*VL*a - 2*c*R0**2*VL + c*R0*VH*a**2 + 4*c*R0*VH*a + 4*c*R0*VH + c*R0*VL*a**2 + 4*c*R0*VL*a + 4*c*R0*VL + R0**2*VH*VL*a + 2*R0**2*VH*VL)) + n - (VH - VL)*(n - (-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)))*(a + 1)*(0.5*x - 0.5*(a + 1)/(a + 2))/(a + 2)**2 + (x - (a + 1)/(a + 2))*(-c*(a + 1)*(0.5*x - 0.5*(a + 1)/(a + 2))/R0 + (0.5*n - 0.5*(-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)))*(a + 2)*(c*(-R0 + a + 2) + R0*VH)*(c*(-R0 + a + 2) + R0*VL)/(R0**2*(VH - VL)**2)) - (-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)),0)
# plt.contour(x, n, -R0*(x - (a + 1)/(a + 2))*(c*(a + 1)*(VH*a + 2*VH - VL*a - 2*VL) - sqrt(-(VH - VL)*(a + 1)*(a + 2)*(4*c**2*R0**2 - 8*c**2*R0*a - 16*c**2*R0 - c**2*VH*a**2 - 3*c**2*VH*a - 2*c**2*VH + c**2*VL*a**2 + 3*c**2*VL*a + 2*c**2*VL + 4*c**2*a**2 + 16*c**2*a + 16*c**2 - 4*c*R0**2*VH - 4*c*R0**2*VL + 4*c*R0*VH*a + 8*c*R0*VH + 4*c*R0*VL*a + 8*c*R0*VL + 4*R0**2*VH*VL)))*(VH**2 - 2*VH*VL + VL**2)/(2*(VH*a + 2*VH - VL*a - 2*VL)*(c**2*R0**2*a + 2*c**2*R0**2 - 2*c**2*R0*a**2 - 8*c**2*R0*a - 8*c**2*R0 + c**2*a**3 + 6*c**2*a**2 + 12*c**2*a + 8*c**2 - c*R0**2*VH*a - 2*c*R0**2*VH - c*R0**2*VL*a - 2*c*R0**2*VL + c*R0*VH*a**2 + 4*c*R0*VH*a + 4*c*R0*VH + c*R0*VL*a**2 + 4*c*R0*VL*a + 4*c*R0*VL + R0**2*VH*VL*a + 2*R0**2*VH*VL)) + n - (-c*(-R0 + a + 2) - R0*VL)/(R0*(VH - VL)),0, colors='orange', linewidths=2, linestyles='dashdot')

font1 = {'family': 'SimHei', 'weight': 'bold', 'size': 30}
plt.annotate(text='unstable limit cycle',xy=(0.8,0.05),xytext=(0.7,0.25), arrowprops=dict(facecolor='black', edgecolor='darkgrey',headlength=15, width=5, headwidth=18), font=font1, horizontalalignment='center', verticalalignment='center')
# plt.annotate(text='stable limit cycle',xy=(0.5,0.91),xytext=(0.5,0.75), arrowprops=dict(facecolor='black', edgecolor='darkgrey',headlength=15, width=5, headwidth=18), font=font1, horizontalalignment='center', verticalalignment='center')

# plt.annotate(text='',xytext=(1 - c/(R0 * (c - VH)) + 0.003, 0.98),xy=(1 - c/(R0 * (c - VH)), 1), arrowprops=dict(facecolor='darkgrey', edgecolor='darkgrey',headlength=15,width=10))
# plt.scatter((a + 1) / (a + 2), -1 * (c * (2 - R0 + a) + R0 * VL) / (R0 * (VH - VL)),color='white', marker='o', edgecolors='black', s=200, linewidths=2)
# plt.scatter(1 - c/(R0 * (c - VL)), 0, c='blue', s=200, edgecolors='black', linewidths=2)
# plt.scatter(1 - c/(R0 * (c - VH)), 1, c='blue', s=200, edgecolors='black', linewidths=2)
plt.yticks(fontproperties='Times New Roman', size=30)
plt.xticks(fontproperties='Times New Roman', size=30)
plt.xlim((-0.02, 1.02))
plt.ylim((-0.02, 1.02))
plt.xlabel('Fraction of cooperators x', font)
plt.ylabel('Environment n', font)
plt.title('General imitation model with Fermi function\n', fontproperties='Times New Roman', size = 35)

# plt.savefig("fig2-1.png")
plt.show()
