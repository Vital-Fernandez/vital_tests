# import matplotlib
# matplotlib.use('QT4Agg')
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import axes3d
#
# class MyAxes3D(axes3d.Axes3D):
#
#     def __init__(self, baseObject, sides_to_draw):
#         self.__class__ = type(baseObject.__class__.__name__,
#                               (self.__class__, baseObject.__class__),
#                               {})
#         self.__dict__ = baseObject.__dict__
#         self.sides_to_draw = list(sides_to_draw)
#         self.mouse_init()
#
#     def set_some_features_visibility(self, visible):
#         for t in self.w_zaxis.get_ticklines() + self.w_zaxis.get_ticklabels():
#             t.set_visible(visible)
#         self.w_zaxis.line.set_visible(visible)
#         self.w_zaxis.pane.set_visible(visible)
#         self.w_zaxis.label.set_visible(visible)
#
#     def draw(self, renderer):
#         # set visibility of some features False
#         self.set_some_features_visibility(False)
#         # draw the axes
#         super(MyAxes3D, self).draw(renderer)
#         # set visibility of some features True.
#         # This could be adapted to set your features to desired visibility,
#         # e.g. storing the previous values and restoring the values
#         self.set_some_features_visibility(True)
#
#         zaxis = self.zaxis
#         draw_grid_old = zaxis.axes._draw_grid
#         # disable draw grid
#         zaxis.axes._draw_grid = False
#
#         tmp_planes = zaxis._PLANES
#
#         if 'l' in self.sides_to_draw :
#             # draw zaxis on the left side
#             zaxis._PLANES = (tmp_planes[2], tmp_planes[3],
#                              tmp_planes[0], tmp_planes[1],
#                              tmp_planes[4], tmp_planes[5])
#             zaxis.draw(renderer)
#         if 'r' in self.sides_to_draw :
#             # draw zaxis on the right side
#             zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
#                              tmp_planes[1], tmp_planes[0],
#                              tmp_planes[4], tmp_planes[5])
#             zaxis.draw(renderer)
#
#         zaxis._PLANES = tmp_planes
#
#         # disable draw grid
#         zaxis.axes._draw_grid = draw_grid_old
#
# def example_surface(ax):
#     """ draw an example surface. code borrowed from http://matplotlib.org/examples/mplot3d/surface3d_demo.html """
#     from matplotlib import cm
#     import numpy as np
#     X = np.arange(-5, 5, 0.25)
#     Y = np.arange(-5, 5, 0.25)
#     X, Y = np.meshgrid(X, Y)
#     R = np.sqrt(X**2 + Y**2)
#     Z = np.sin(R)
#     surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
#
# if __name__ == '__main__':
#     fig = plt.figure(figsize=(15, 5))
#     ax = fig.add_subplot(131, projection='3d')
#     ax.set_title('z-axis left side')
#     ax = fig.add_axes(MyAxes3D(ax, 'l'))
#     example_surface(ax) # draw an example surface
#     ax = fig.add_subplot(132, projection='3d')
#     ax.set_title('z-axis both sides')
#     ax = fig.add_axes(MyAxes3D(ax, 'lr'))
#     example_surface(ax) # draw an example surface
#     ax = fig.add_subplot(133, projection='3d')
#     ax.set_title('z-axis right side')
#     ax = fig.add_axes(MyAxes3D(ax, 'r'))
#     # example_surface(ax) # draw an example surface
#     plt.show()

# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import src.specsiser as sr
from pathlib import Path
from matplotlib import pyplot as plt, rcParams, colors, cm, ticker, gridspec
from matplotlib.cbook import get_sample_data
from matplotlib.offsetbox   import OffsetImage, AnnotationBbox
from astro.data.muse.common_methods import STANDARD_AXES, DARK_PLOT, background_color, foreground_color

from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import numpy as np

# %matplotlib inline

defaultConf = DARK_PLOT.copy()
rcParams.update(defaultConf)


plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['font.sans-serif'] = ['SimHei']

def plot_opaque_cube(x=10, y=20, z=30, dx=40, dy=50, dz=60):

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')


    xx = np.linspace(x, x+dx, 2)
    yy = np.linspace(y, y+dy, 2)
    zz = np.linspace(z, z+dz, 2)

    xx, yy = np.meshgrid(xx, yy)

    ax.plot_surface(xx, yy, z)
    ax.plot_surface(xx, yy, z+dz)

    yy, zz = np.meshgrid(yy, zz)
    ax.plot_surface(x, yy, zz)
    ax.plot_surface(x+dx, yy, zz)

    xx, zz = np.meshgrid(xx, zz)
    ax.plot_surface(xx, y, zz)
    ax.plot_surface(xx, y+dy, zz)
    # ax.set_xlim3d(-dx, dx*2, 20)
    # ax.set_xlim3d(-dx, dx*2, 20)
    # ax.set_xlim3d(-dx, dx*2, 20)
    plt.title("Cube")
    plt.show()


def plot_linear_cube(x, y, z, dx, dy, dz, color='red'):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.set_xlabel('$log(U)$')
    ax.set_ylabel('$T_{eff}$')
    ax.set_zlabel('$12+log(O/H)$')
    xx = [x, x, x+dx, x+dx, x]
    yy = [y, y+dy, y+dy, y, y]
    kwargs = {'alpha': 0, 'color': color}
    ax.plot3D(xx, yy, [z]*5, **kwargs)
    ax.plot3D(xx, yy, [z+dz]*5, **kwargs)
    ax.plot3D([x, x], [y, y], [z, z+dz], **kwargs)
    ax.plot3D([x, x], [y+dy, y+dy], [z, z+dz], **kwargs)
    ax.plot3D([x+dx, x+dx], [y+dy, y+dy], [z, z+dz], **kwargs)
    ax.plot3D([x+dx, x+dx], [y, y], [z, z+dz], **kwargs)
    plt.title('Cube')
    plt.show()


if __name__ == "__main__":
    plot_linear_cube(30000, -4.0, 7.1, 90000, -1.5, 8.9)
