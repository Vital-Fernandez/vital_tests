import matplotlib.pyplot as plt
import numpy as np
from matplotlib.widgets  import RectangleSelector

xdata = np.linspace(0,9*np.pi, num=301)
ydata = np.sin(xdata)

fig, ax = plt.subplots()
line, = ax.plot(xdata, ydata)


def line_select_callback(eclick, erelease):
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata

    rect = plt.Rectangle( (min(x1,x2),min(y1,y2)), np.abs(x1-x2), np.abs(y1-y2) )
    ax.add_patch(rect)


rs = RectangleSelector(ax, line_select_callback,
                       drawtype='box', useblit=False, button=[1],
                       minspanx=5, minspany=5, spancoords='pixels',
                       interactive=True)

plt.show()

# import logging
# import matplotlib
# from matplotlib.widgets import Lasso
# from matplotlib.colors import colorConverter
# from matplotlib.collections import RegularPolyCollection
# from matplotlib import path
# import numpy as np
# import matplotlib.pyplot as plt
# from numpy.random import rand
#
# logger = logging.getLogger()
# logger.setLevel(logging.DEBUG)
#
#
# class Datum(object):
#       colorin = colorConverter.to_rgba('red')
#       colorShift = colorConverter.to_rgba('cyan')
#       colorCtrl = colorConverter.to_rgba('pink')
#       colorout = colorConverter.to_rgba('blue')
#
#       def __init__(self, x, y, include=False):
#          self.x = x
#          self.y = y
#          if include:
#             self.color = self.colorin
#          else:
#             self.color = self.colorout
#
#
# class LassoManager(object):
#     def __init__(self, ax, data):
#         self.axes = ax
#         self.canvas = ax.figure.canvas
#         self.data = data
#
#         self.Nxy = len(data)
#
#         facecolors = [d.color for d in data]
#         self.xys = [(d.x, d.y) for d in data]
#         fig = ax.figure
#         self.collection = RegularPolyCollection(fig.dpi, 6, sizes=(100,), facecolors=facecolors, offsets=self.xys,
#                                                 transOffset=ax.transData)
#
#         ax.add_collection(self.collection)
#
#         self.cid = self.canvas.mpl_connect('button_press_event', self.onpress)
#         self.keyPress = self.canvas.mpl_connect('key_press_event', self.onKeyPress)
#         self.keyRelease = self.canvas.mpl_connect('key_release_event', self.onKeyRelease)
#         self.pick = self.canvas.mpl_connect('pick_event', self.onpick)
#         self.lasso = None
#         self.shiftKey = False
#         self.ctrlKey = False
#         self.pickEvent = False
#
#     def callback(self, verts):
#         logging.debug('in LassoManager.callback(). Shift: %s, Ctrl: %s' % (self.shiftKey, self.ctrlKey))
#         facecolors = self.collection.get_facecolors()
#         p = path.Path(verts)
#         ind = p.contains_points(self.xys)
#         for i in range(len(self.xys)):
#             if ind[i]:
#                 if self.shiftKey:
#                     facecolors[i] = Datum.colorShift
#                 elif self.ctrlKey:
#                     facecolors[i] = Datum.colorCtrl
#                 else:
#                     facecolors[i] = Datum.colorin
#                     print
#                     self.xys[i]
#             else:
#                 facecolors[i] = Datum.colorout
#
#         self.canvas.draw_idle()
#         self.canvas.widgetlock.release(self.lasso)
#         del self.lasso
#
#     def onpress(self, event):
#         logging.debug('in LassoManager.onpress(). Event received: %s' % event)
#         if self.pickEvent:
#             self.pickEvent = False
#             return
#         if self.canvas.widgetlock.locked(): return
#         if event.inaxes is None: return
#         self.lasso = Lasso(event.inaxes, (event.xdata, event.ydata), self.callback)
#         # acquire a lock on the widget drawing
#         self.canvas.widgetlock(self.lasso)
#
#     def onKeyPress(self, event):
#         logging.debug('in LassoManager.onKeyPress(). Event received: %s (key: %s)' % (event, event.key))
#         if event.key == 'alt':
#             self.ctrlKey = True
#         if event.key == 'shift':
#             self.shiftKey = True
#
#     def onKeyRelease(self, event):
#         logging.debug('in LassoManager.onKeyRelease(). Event received: %s (key: %s)' % (event, event.key))
#         if event.key == 'alt':
#             self.ctrlKey = False
#         if event.key == 'shift':
#             self.shiftKey = False
#
#     def onpick(self, event):
#         logging.debug('in LassoManager.onpick(). Event received: %s' % event)
#         self.pickEvent = True
#         if event.mouseevent.button == 3:
#             index = event.ind
#             print('onpick scatter: ', index, np.take(x, index), np.take(y, index))
#
#
# if __name__ == '__main__':
#     x,y =rand(2,100)
#     data = [Datum(*xy) for xy in zip(x,y)]
#     fig = plt.figure()
#     ax = plt.axes()
#
#     ax.scatter(x,y,picker=True)
#
#     lman = LassoManager(ax, data)
#     plt.show()