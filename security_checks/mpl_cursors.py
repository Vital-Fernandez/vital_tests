import numpy as np
import matplotlib.pyplot as plt
import mplcursors

x_1 = np.array([1, 2, 3])
y_1 = x_1
labels1 = ['a', 'b', 'c']

x_2 = np.array([1, 2, 3, 4])
y_2 = x_2 * 2
labels2 = ['d', 'e', 'f', 'g']

x_values = [x_1, x_2]
y_values = [y_1, y_2]
labels = [labels1, labels2]

fig, ax = plt.subplots(figsize=(10, 10))

for i, x in enumerate(x_values):
    data_scatter = ax.scatter(x, y_values[i], label=f'Data set {i+1}')
    labels_group = labels[i]
    mplcursors.cursor(data_scatter).connect("add", lambda sel, labels_group=labels_group: sel.annotation.set_text(labels_group[sel.target.index]))

ax.legend()
plt.show()
