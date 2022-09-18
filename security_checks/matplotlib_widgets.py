import matplotlib.pyplot as plt
plt.rcParams['toolbar'] = 'toolmanager'
from matplotlib.backend_tools import ToolBase, ToolToggleBase


class ListTools(ToolBase):
    """List all the tools controlled by the `ToolManager`."""
    # keyboard shortcut
    default_keymap = 'm'
    description = 'List Tools'

    def trigger(self, *args, **kwargs):
        print('hola')
        print('_' * 80)
        print("bicho{0:12} {1:45} {2}".format(
            'Name (id)', 'Tool description', 'Keymap'))
        print('-' * 80)
        tools = self.toolmanager.tools
        for name in sorted(tools):
            if not tools[name].description:
                continue
            keys = ', '.join(sorted(self.toolmanager.get_tool_keymap(name)))
            print("bicho{0:12} {1:45} {2}".format(
                name, tools[name].description, keys))
        print('_' * 80)
        print("Active Toggle tools")
        print("{0:12} {1:45}".format("Group", "Active"))
        print('-' * 80)
        for group, active in self.toolmanager.active_toggle.items():
            print("{0:12} {1:45}".format(str(group), str(active)))

        for i in range(10):
            print(i)


class GroupHideTool(ToolToggleBase):
    """Show lines with a given gid."""
    default_keymap = 'S'
    description = 'Show by gid'
    default_toggled = True

    def __init__(self, *args, gid, **kwargs):
        self.gid = gid
        super().__init__(*args, **kwargs)

    def enable(self, *args):
        self.set_lines_visibility(True)

    def disable(self, *args):
        self.set_lines_visibility(False)

    def set_lines_visibility(self, state):
        for ax in self.figure.get_axes():
            for line in ax.get_lines():
                if line.get_gid() == self.gid:
                    line.set_visible(state)
        self.figure.canvas.draw()


fig = plt.figure()
plt.plot([1, 2, 3], gid='mygroup')
plt.plot([2, 3, 4], gid='unknown')
plt.plot([3, 2, 1], gid='mygroup')

# Add the custom tools that we created
fig.canvas.manager.toolmanager.add_tool('List', ListTools)
fig.canvas.manager.toolmanager.add_tool('Show', GroupHideTool, gid='mygroup')

fig.canvas.manager.toolbar.add_tool('Show', 'navigation', 1)

plt.show()

# import matplotlib.pyplot as plt
# from matplotlib.backend_tools import Cursors
#
#
# fig, axs = plt.subplots(len(Cursors), figsize=(6, len(Cursors) + 0.5),
#                         gridspec_kw={'hspace': 0})
# fig.suptitle('Hover over an Axes to see alternate Cursors')
#
# for cursor, ax in zip(Cursors, axs):
#     ax.cursor_to_use = cursor
#     ax.text(0.5, 0.5, cursor.name,
#             horizontalalignment='center', verticalalignment='center')
#     ax.set(xticks=[], yticks=[])
#
#
# def hover(event):
#     if fig.canvas.widgetlock.locked():
#         # Don't do anything if the zoom/pan tools have been enabled.
#         return
#
#     fig.canvas.set_cursor(
#         event.inaxes.cursor_to_use if event.inaxes else Cursors.POINTER)
#
#
# fig.canvas.mpl_connect('motion_notify_event', hover)
#
# plt.show()
# import matplotlib.pyplot as plt
# from matplotlib.backend_bases import NavigationToolbar2, Event
#
# home = NavigationToolbar2.home
#
# def new_home(self, *args, **kwargs):
#     s = 'home_event'
#     event = Event(s, self)
#     event.foo = 100
#     self.canvas.callbacks.process(s, event)
#     home(self, *args, **kwargs)
#
# NavigationToolbar2.home = new_home
#
# def handle_home(evt):
#     print('new home')
#     print(evt.foo)
#
# fig = plt.figure()
# fig.canvas.mpl_connect('home_event', handle_home)
# plt.text(0.35, 0.5, 'Hello world!', dict(size=30))
# plt.show()


# import matplotlib.pyplot as plt
# import numpy as np
# from matplotlib.widgets  import RectangleSelector
#
# xdata = np.linspace(0,9*np.pi, num=301)
# ydata = np.sin(xdata)
#
# fig, ax = plt.subplots()
# line, = ax.plot(xdata, ydata)
#
#
# def line_select_callback(eclick, erelease):
#     x1, y1 = eclick.xdata, eclick.ydata
#     x2, y2 = erelease.xdata, erelease.ydata
#
#     rect = plt.Rectangle( (min(x1,x2),min(y1,y2)), np.abs(x1-x2), np.abs(y1-y2) )
#     ax.add_patch(rect)
#
#
# rs = RectangleSelector(ax, line_select_callback,
#                        drawtype='box', useblit=False, button=[1],
#                        minspanx=5, minspany=5, spancoords='pixels',
#                        interactive=True)
#
# plt.show()

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