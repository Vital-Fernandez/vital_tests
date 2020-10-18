from matplotlib import pyplot as plt, gridspec

fig = plt.figure()
gs = gridspec.GridSpec(nrows=1,
                       ncols=2,
                       figure=fig,
                       width_ratios=[1, 2],
                       height_ratios=[1])
ax1 = fig.add_subplot(gs[0])
ax1.text(0.5, 0.5, 'ax1: gs[0]', fontsize=12, fontweight="bold", va="center", ha="center")
ax2 = fig.add_subplot(gs[1])
ax2.text(0.5, 0.5, 'ax2: gs[1]', fontsize=12, fontweight="bold", va="center", ha="center")
plt.show()

# fig = plt.figure(figsize=(7,7))
# gs = gridspec.GridSpec(nrows=3,
#                        ncols=3,
#                        figure=fig,
#                        width_ratios= [1, 1, 1],
#                        height_ratios=[1, 1, 1],
#                        wspace=0.3,
#                        hspace=0.3)
# ax1 = fig.add_subplot(gs[0, 0])
# ax1.text(0.5, 0.5, 'ax1: gs[0, 0]', fontsize=12, fontweight="bold", va="center", ha="center")  # adding text to ax1
# ax2 = fig.add_subplot(gs[0, 1:3])
# ax2.text(0.5, 0.5, 'ax2: gs[0, 1:3]', fontsize=12, fontweight="bold", va="center", ha="center")
# ax3 = fig.add_subplot(gs[1:3, 0:2])
# ax3.text(0.5, 0.5, 'ax3: gs[1:3, 0:2]', fontsize=12, fontweight="bold", va="center", ha="center")
# ax4 = fig.add_subplot(gs[1:3, 2])
# ax4.text(0.5, 0.5, 'ax4: gs[1:3, 2]', fontsize=12, fontweight="bold", va="center", ha="center")
# plt.show()