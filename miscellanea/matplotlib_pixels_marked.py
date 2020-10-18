import numpy as np
import matplotlib.pyplot as plt

data = np.arange(100).reshape(10, 10)

cmap = plt.cm.gray
norm = plt.Normalize(data.min(), data.max())
rgba = cmap(norm(data))

# Set the diagonal to red
rgba[range(10), range(10), :3] = 1, 0, 0

fig, ax = plt.subplots()
ax.imshow(rgba, interpolation='nearest')

# Add the colorbar using a fake (not shown) image.
im = ax.imshow(data, visible=False, cmap=cmap)
fig.colorbar(im)

plt.show()