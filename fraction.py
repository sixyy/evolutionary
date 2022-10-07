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
from scipy import integrate
import sympy as sm
import pandas as pd    #Processing Data from Spreadsheet
from scipy.integrate import odeint #Ordinary differential equation solver
import matplotlib.pyplot as plt

# x, n ,p= sm.symbols('x n p')
x= sm.symbols('x')
n = sm.symbols('n')
p= sm.symbols('p')
# fig4-3-1
R0 = 6.0
S0 = 6.0
T0 = 4.0
P0 = 1.0
R1 = 2.0
S1 = 4.0
T1 = 3.0
P1 = 6.0
theta = 1.0
episi = 1.0
beta = 1.0
y0 = [0.9, 0.11]
y1 = [0.2, 0.04]

# fig4-3-2
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
# y1 = [0.5, 0.2]
t = np.linspace(0, 1000, 10000)

def pend(y, t, R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta):
    x, n = y
    dydt = [x*(1 - x)*(1/(1 + math.exp(-beta* ((((R1 - R0 - S1 + S0)*n + (R0 - S0))*x + (S1 - S0)*n + S0) - (((T1 - T0 - P1 + P0)*n + (T0 - P0))* x + (P1 - P0)*n + P0)))) - 1/(1 + math.exp(beta *((((R1 - R0 - S1 + S0)*n + (R0 - S0))*x + (S1 - S0)*n + S0) - (((T1 - T0 - P1 + P0)*n + (T0 - P0))*x + (P1 - P0)*n + P0))))), episi*n*(1 - n)*(-1 + (1 + theta)*x)]
    return dydt

sol = odeint(pend, y0, t, args=(R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta))
sol1 = odeint(pend, y1, t, args=(R0, S0, T0, P0, R1, S1, T1, P1, theta, episi, beta))
# print(t,sol[:,0])
df = pd.DataFrame({
    'day': t,
    'x with initial point (0.9,0.11)': sol[:,0],
    'n with initial point (0.9,0.11)': sol[:,1],
    'x with initial point (0.2,0.04)': sol1[:,0],
    'n with initial point (0.2,0.04)': sol1[:,1]
})

font = {'family': 'Times New Roman', 'weight': 'normal', 'size': 30}

# ':', '--', '-', '-.', '-'
# , '#cc8ac0', 'grey', 'green'
#plt.style.use('ggplot')
df.plot(x='day',
        y=['x with initial point (0.2,0.04)', 'n with initial point (0.2,0.04)', 'x with initial point (0.9,0.11)', 'n with initial point (0.9,0.11)'],
        color=['#bb6424', 'lightgreen', 'grey', '#aac6ca'],
        style=[':', '--', '-', '-'],
        linewidth=3,
        figsize=(10, 10),
        fontsize=15)
# print(x)
# df.loc['Col_sum'] = df.apply(lambda x: x.sum())
# m = df.at['Col_sum','Vaccinated']
# plt.axhline(m/N, linewidth=5, color='grey')
# plt.annotate(text='average vaccination level',xytext=(60, m/N + 0.3),xy=(60, m/N), arrowprops=dict(facecolor='black', edgecolor='darkgrey',headlength=15, width=5, headwidth=18), font=font, horizontalalignment='center', verticalalignment='center')

# df.drop('Col_sum', inplace = True)
# plt.text(15, 0.5, r'\theta=0.01', font, wrap=True)

ax = plt.gca() #获取边框
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)



plt.xlim((0, 1000))
plt.ylim((-0.02, 1.02))
# plt.grid(ls=':')


plt.xlabel('Time', font)
plt.ylabel('System state', font)
plt.legend(fontsize=14.5, loc='upper right')
# plt.tight_layout()
# print(x)
# plt.axhline(0.64446901, linewidth=5, color='grey')

plt.savefig("fig5.pdf")

plt.show()