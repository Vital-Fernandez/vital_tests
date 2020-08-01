import numpy as np
import matplotlib.pyplot as plt


def onclick(event):
    print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          ('double' if event.dblclick else 'single', event.button, event.x, event.y, event.xdata, event.ydata))


fig, ax = plt.subplots()
ax.plot(np.random.rand(10))
cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.show()


def Span_Manager(Wlow, Whig):
    # Check to make sure we are at the measuring lines screen
    if (pv.CurrentScreen == 'LINEMESURER'):

        # Check we are not just clicking on the plot
        if (Wlow != Whig):

            # Clear the plot
            pv.reset_fig()

            # Case selecting 1/3 region
            if len(pv.Selections) == 0:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)

            # Case selecting 2/3 region
            elif len(pv.Selections) == 2:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)
                pv.Selections.sort()

            # Case selecting 3/3 region
            elif len(pv.Selections) == 4:
                pv.Selections.append(Wlow)
                pv.Selections.append(Whig)
                pv.Selections.sort()

            elif len(pv.Selections) == 6:

                # Caso que se corrija la region de la linea
                if (Wlow > pv.Selections[1] and Whig < pv.Selections[4]):
                    pv.Selections[2] = Wlow
                    pv.Selections[3] = Whig

                # Caso que se corrija el continuum izquierdo
                elif (Wlow < pv.Selections[2] and Whig < pv.Selections[2]):
                    pv.Selections[0] = Wlow
                    pv.Selections[1] = Whig

                # Caso que se corrija el continuum derecho
                elif (Wlow > pv.Selections[3] and Whig > pv.Selections[3]):
                    pv.Selections[4] = Wlow
                    pv.Selections[5] = Whig

                # Case we want to select the complete region
                elif (Wlow < pv.Selections[0] and Whig > pv.Selections[5]):

                    # Remove line from dataframe and save it
                    pv.remove_lines_df(pv.current_df, pv.Current_Label)

                    # Save lines log df
                    pv.save_lineslog_dataframe(pv.current_df, pv.lineslog_df_address)

                    # Clear the selections
                    del pv.Selections[:]

            elif len(pv.Selections) > 6:
                print
                "WARNING: Selections size over limits 6"

            pv.ManagePlots()

        pv.FigCanvas.draw()

    return