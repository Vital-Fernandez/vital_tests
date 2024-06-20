import matplotlib.pyplot as plt
from matplotlib import rc_context
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import lime

lime.theme.set_style('dark')
fig_cfg = lime.theme.fig_defaults()

fig_cfg['xtick.bottom'] = False
fig_cfg['xtick.labelbottom'] = False
fig_cfg['xtick.top'] = False
fig_cfg['xtick.labeltop'] = False

fig_cfg['ytick.left'] = False
fig_cfg['ytick.labelleft'] = False
fig_cfg['ytick.right'] = False
fig_cfg['ytick.labelright'] = False





# Create axis
axes = [5, 5, 5]

# Create Data
data = np.ones(axes, dtype=np.bool_)

# Control Transparency
alpha = 0.9

# Control colour
colors = np.empty(axes + [4], dtype=np.float32)

colors = 'red'

with rc_context(fig_cfg):
    # Plot figure
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Voxels is used to customizations of the
    # sizes, positions and colors.
    ax.voxels(data, facecolors=colors)
    plt.tight_layout()
    plt.savefig('/home/vital/Desktop/red_cube.png', transparent=True)
