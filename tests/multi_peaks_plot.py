from matplotlib import pyplot as plt, rc_context
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from lime import theme

theme.set_style('dark', fig_cfg={"figure.dpi": 300})

fig_cfg = theme.fig_defaults()

with rc_context(fig_cfg):

    # Creating a figure with multiple peaks and a transparent background
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    X = np.linspace(-5, 5, 100)
    Y = np.linspace(-5, 5, 100)
    X, Y = np.meshgrid(X, Y)
    Z = np.sin(np.sqrt(X**2 + Y**2)) + np.cos(X) + np.sin(Y)

    surf = ax.plot_surface(X, Y, Z, cmap='coolwarm')

    # Customizing the font color to white for axis ticks and labels
    # ax.xaxis.set_tick_params(labelcolor='white')
    # ax.yaxis.set_tick_params(labelcolor='white')
    # ax.zaxis.set_tick_params(labelcolor='white')
    #
    # ax.set_xlabel('X axis', color='white')
    # ax.set_ylabel('Y axis', color='white')
    # ax.set_zlabel('Z axis', color='white')

    # Removing tick labels
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    # Setting the background to transparent
    fig.patch.set_alpha(0)

    plt.show()