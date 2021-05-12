import numpy as np
import matplotlib.pyplot as plt

import matplotlib.pyplot as plt

plt.plot(range(10))

def onclick(event):
    print(event.key)

cid = plt.gcf().canvas.mpl_connect('key_press_event', onclick)

plt.show()

# class push_to_advance(object):
#
#     def __init__(self, input_list):
#         self.fig = plt.figure()
#         self.ax = self.fig.gca()
#         self.bound_keys = []
#         self.bound_cid = {}
#         self.objList = input_list
#
#     def add_step_through(self, gen, key):
#         key = key[0] # make a single char
#         if key in self.bound_keys:
#             raise RuntimeError("key %s already bound"%key)
#         idx_data, first_data = next(gen)
#         self.objList[idx_data] = key
#         self.ax.plot(first_data)
#         self.ax.update({'title': idx_data, 'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'$\gamma_{\nu} (erg\,cm^{3} s^{-1})$'})
#         self.fig.canvas.draw()
#         self.bound_keys.append(key)
#
#         def ontype(event):
#             if event.key == key:
#                 try:
#                     plt.cla()
#                     idx_data, secondary_data = next(gen)
#                     self.objList[idx_data] = key
#                     self.ax.plot(secondary_data)
#                     self.ax.update({'title': idx_data, 'xlabel': r'Wavelength $(\AA)$', 'ylabel': r'$\gamma_{\nu} (erg\,cm^{3} s^{-1})$'})
#                     self.fig.canvas.draw()
#                 except StopIteration:
#                     self.fig.canvas.mpl_disconnect(self.bound_cid[key])
#                     del self.bound_cid[key]
#                     self.bound_keys.remove(key)
#
#         self.bound_cid[key] = self.fig.canvas.mpl_connect('key_press_event', ontype)
#
#
# # Define the event
# def ugly_math():
#     print('you will hit this once')
#     for j in range(10):
#         print('loop ', j)
#         # insert math here
#         yield j, np.random.random(10) * j
#
# fun = ugly_math()
# my_list = [None] * 15
# print(my_list)
#
# pta = push_to_advance(my_list)
# gen = ugly_math()
# pta.add_step_through(gen, 'a')
# plt.show()
# print(my_list)
# print(pta.objList)
# # test_array = np.arange(100).reshape(10, 10)
# # pta.add_step_through(test_array.__iter__(), 'b')