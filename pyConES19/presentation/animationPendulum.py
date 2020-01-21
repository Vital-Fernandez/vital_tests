import math
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
import matplotlib.pyplot as plt
import os

background = np.array((43, 43, 43))/255.0
foreground = np.array((179, 199, 216))/255.0
red = np.array((43, 43, 43))/255.0
yellow = np.array((191, 144, 0))/255.0

figConf = {'text.color': foreground,
            'figure.figsize': (10,10),
            'figure.facecolor':background,
            'axes.facecolor':background,
            'axes.edgecolor':foreground,
            'axes.labelcolor':foreground,
            'axes.labelsize':18,
            'xtick.labelsize':16,
            'ytick.labelsize':16,
            'xtick.color':foreground,
            'ytick.color':foreground,
            'legend.edgecolor':'inherit',
            'legend.facecolor':'inherit',
            'legend.fontsize':16,
             'legend.loc':"center right"}

matplotlib.rcParams.update(figConf)

# Initial and end values
st = 0          # Start time (s)
et = 20.4       # End time (s)
ts = 0.1        # Time step (s)
g = 9.81        # Acceleration due to gravity (m/s^2)
L = 1           # Length of pendulum (m)
b = 0.5         # Damping factor (kg/s)
m = 1           # Mass of bob (kg)

dataFolder = 'E:/Dropbox/Astrophysics/Seminars/PyConES_2019/pendulum/'

# 1st order equations to solve in a function
"""
 theta1 is angular displacement at current time instant
 theta2 is angular velocity at current time instant
 dtheta2_dt is angular acceleration at current time instant
 dtheta1_dt is rate of change of angular displacement at current time instant i.e. same as theta2 
"""

def sim_pen_eq(t, theta):
    dtheta2_dt = (-b/m)*theta[1] + (-g/L)*np.sin(theta[0])
    dtheta1_dt = theta[1]
    return [dtheta1_dt, dtheta2_dt]

# main
theta1_ini = 0                 # Initial angular displacement (rad)
theta2_ini = 3                 # Initial angular velocity (rad/s)
theta_ini = [theta1_ini, theta2_ini]
t_span = [st,et+ts]
t = np.arange(st,et+ts,ts)
sim_points = len(t)
l = np.arange(0,sim_points,1)

theta12 = solve_ivp(sim_pen_eq, t_span, theta_ini, t_eval = t)
theta1 = theta12.y[0,:]
theta2 = theta12.y[1,:]
# plt.plot(t, theta1, label='Angular Displacement (rad)')
# plt.plot(t, theta2, label='Angular velocity (rad/s)')
# plt.xlabel('Time(s)')
# plt.ylabel('Angular Disp.(rad) and Angular Vel.(rad/s)')
# plt.legend()
# plt.show()


# Simulation
x = L*np.sin(theta1)
y = -L*np.cos(theta1)

for point in l:
    plt.figure()
    plt.plot(x[point],y[point],'o', color=foreground, markersize=20)
    plt.plot([0,x[point]], [0,y[point]], color=foreground)
    plt.xlim(-L-0.5,L+0.5)
    plt.ylim(-L-0.5,L+0.5)
    plt.xlabel('x direction')
    plt.ylabel('y direction')
    plt.tick_params(axis='x', bottom=False, top=False, labelbottom=False)
    plt.tick_params(axis='y', left=False, labelleft=False)
    filenumber = point
    filenumber = format(filenumber, "05")
    filename = dataFolder + "image{}.png".format(filenumber)
    plt.savefig(filename, bbox_inches='tight', facecolor=background)
    plt.close()


# import numpy as np
# import matplotlib
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation
# from scipy.integrate import odeint
#
# background = np.array((43, 43, 43))/255.0
# foreground = np.array((179, 199, 216))/255.0
# red = np.array((43, 43, 43))/255.0
# yellow = np.array((191, 144, 0))/255.0
#
# figConf = {'text.color': foreground,
#             'figure.figsize': (10,10),
#             'figure.facecolor':background,
#             'axes.facecolor':background,
#             'axes.edgecolor':foreground,
#             'axes.labelcolor':foreground,
#             'axes.labelsize':18,
#             'xtick.labelsize':16,
#             'ytick.labelsize':16,
#             'xtick.color':foreground,
#             'ytick.color':foreground,
#             'legend.edgecolor':'inherit',
#             'legend.facecolor':'inherit',
#             'legend.fontsize':16,
#              'legend.loc':"center right"}
#
# matplotlib.rcParams.update(figConf)
#
#
# def pend(y, t, b, c):
#     theta, omega = y
#     dydt = [omega, -b*omega - c*np.sin(theta)]
#     return dydt
#
# b = 0.25
# c = 5.0
# y0 = [np.pi - 0.1, 0.0]
# t = np.linspace(0, 20, 301)
#
# sol = odeint(pend, y0, t, args=(b, c))
#
# fig = plt.figure(figsize=(5, 5), facecolor='w')
# ax = fig.add_subplot(1, 1, 1)
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.tick_params(which='both', axis='x', bottom=False, top=False, labelbottom=False)
# ax.tick_params(which='both', axis='y', left=False, labelleft=False)
#
#
# lns = []
# for i in range(len(sol)):
#     ln, = ax.plot([0, np.sin(sol[i, 0])], [0, -np.cos(sol[i, 0])],
#                   color=foreground, lw=2)
#     # tm = ax.text(-1, 0.9, 'time = %.1fs' % t[i])
#     lns.append([ln])
# ax.set_aspect('equal', 'datalim')
# ani = animation.ArtistAnimation(fig, lns, interval=50)
#
# plt.show()
#
# # ani.save(fn+'.mp4',writer='ffmpeg',fps=1000/50)
# # ani.save(fn+'.gif',writer='imagemagick',fps=1000/50)