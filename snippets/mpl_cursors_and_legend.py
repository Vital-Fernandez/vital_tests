import numpy as np
import matplotlib.pyplot as plt
import mplcursors

x_1 = np.array([1, 2, 3])
y_1 = x_1

x_2 = np.array([1, 2, 3, 4])
y_2 = x_2 * 2

x_values = [x_1, x_2]
y_values = [y_1, y_2]
labels_mplcursors = ['m = 1', 'm = 2']
labels_legend = ['Dataset 1', 'Dataset 2']

fig, ax = plt.subplots(figsize=(10, 10))

container = {}
for i, x in enumerate(x_values):
    line = ax.plot(x, y_values[i], label=labels_legend[i])
    container[labels_mplcursors[i]] = line

for label, line in container.items():
    mplcursors.cursor(line).connect("add", lambda sel, label=label: sel.annotation.set_text(label))


ax.legend()
plt.show()
