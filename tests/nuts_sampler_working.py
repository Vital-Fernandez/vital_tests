import numpy as np
import pymc as pm
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Define a simple Bayesian model
np.random.seed(42)
observed_data = np.random.normal(loc=3.0, scale=1.0, size=100)

with pm.Model() as model:
    mu = pm.Normal('mu', mu=0, sigma=10)
    sigma = pm.HalfNormal('sigma', sigma=1)
    y = pm.Normal('y', mu=mu, sigma=sigma, observed=observed_data)

    # Use NUTS sampler
    trace = pm.sample(1000, tune=500, return_inferencedata=False)

# Create an animation of the sampling process
fig, ax = plt.subplots()

x = np.linspace(-5, 10, 100)
y = np.linspace(0, 0.5, 100)
X, Y = np.meshgrid(x, y)
Z = np.exp(-0.5 * ((X - 3.0) ** 2 / 1.0 ** 2 + Y ** 2))

contour = ax.contourf(X, Y, Z, cmap='viridis')

line, = ax.plot([], [], 'ro', markersize=2)
title = ax.set_title('NUTS Sampling')


def init():
    line.set_data([], [])
    return line, title


def update(frame):
    line.set_data(trace['mu'][:frame], trace['sigma'][:frame])
    title.set_text(f'NUTS Sampling - Iteration {frame}')
    return line, title


ani = FuncAnimation(fig, update, frames=range(1, len(trace['mu'])), init_func=init, blit=True, repeat=False)

plt.show()
