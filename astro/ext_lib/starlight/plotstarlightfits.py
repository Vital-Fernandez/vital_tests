'''
Created on May 11, 2012

@author: william
'''

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astro.ext_lib.starlight.smooth import smooth

class plot_fits_and_SFH(object):
    '''
        This class plots the world most famous graphics of the fits and SFH of Starlight.
    '''

    def __init__(self, sl_out, interactive=True, i_fig=1, figsize=(10, 7)):
        '''
        Constructor
        '''
        self.sl_out = sl_out
        self.fig = plt.figure(i_fig, figsize=figsize)
        if (interactive):
            matplotlib.pyplot.interactive(True)
        else:
            matplotlib.pyplot.interactive(False)

    def set_title(self, title, color='red'):
        self.fig.text(.6, .92, title, color=color)

    def draw_legends(self, metallicity=True, metallicityLoc=2):
        if (metallicity == True): self.axright_upp.legend(loc=metallicityLoc)

    def save(self, outfile='Starlight_fits_and_SFH.png', fmt='png'):
        plt.savefig(outfile, format=fmt)

    def plot_fig_starlight(self, smooth_plot=False):

        self.fig.clf()

        # Some configurations
        Zcolor = ['magenta', 'cyan', 'blue', 'green', 'black', 'red', 'magenta', 'cyan', 'blue', 'green', 'black,',
                  'red']
        lim_res_low = -.29
        lim_res_upp = .99

        # Divide in two boxes 60% and 30% of size.
        box_left_low = [.1, .1, .4, .2]  # left, bottom, width, height
        box_left_upp = [.1, .3, .4, .6]
        box_right_low = [.6, .1, .3, .3]
        box_right_upp = [.6, .4, .3, .3]

        # Create the axes
        self.axleft_low = self.fig.add_axes(box_left_low)
        self.axleft_upp = self.fig.add_axes(box_left_upp, sharex=self.axleft_low)
        self.axright_low = self.fig.add_axes(box_right_low)
        self.axright_upp = self.fig.add_axes(box_right_upp, sharex=self.axright_low)

        # .... Plots ...
        # LEFT
        ### LEFT - UPPER
        if smooth_plot:
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, smooth(self.sl_out.spectra.f_obs, smooth_plot),
                                 color='blue')  # Observed
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, smooth(self.sl_out.spectra.f_syn, smooth_plot),
                                 color='red')  # Synthetic
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, self.sl_out.spectra.f_obs, color='cyan',
                                 alpha=.3)  # Observed
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, self.sl_out.spectra.f_syn, color='magenta',
                                 alpha=.3)  # Synthetic
        else:
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, self.sl_out.spectra.f_obs, color='blue')  # Observed
            self.axleft_upp.plot(self.sl_out.spectra.l_obs, self.sl_out.spectra.f_syn, color='red')  # Synthetic

        np.seterr(divide='ignore')  # Ignore zero-division error.
        err = np.ma.array(np.divide(1., (self.sl_out.spectra.f_wei)))

        if (smooth_plot):
            self.axleft_upp.plot(self.sl_out.spectra.l_obs[err > 0], smooth(err[err > 0], smooth_plot),
                                 color='black')  # Error
            self.axleft_upp.plot(self.sl_out.spectra.l_obs[err > 0], err[err > 0], color='gray', alpha=.3)  # Error
        else:
            self.axleft_upp.plot(self.sl_out.spectra.l_obs[err > 0], err[err > 0], color='black')  # Error

        mi_ = np.min(self.sl_out.spectra.l_obs)
        ma_ = np.max(self.sl_out.spectra.l_obs)
        self.axleft_upp.set_xlim(mi_ * 0.9, ma_ * 1.1)

        ### LEFT - LOWER
        if smooth_plot:
            res = np.ma.array(
                smooth(self.sl_out.spectra.f_obs, smooth_plot) - smooth(self.sl_out.spectra.f_syn, smooth_plot))
        else:
            res = np.ma.array(self.sl_out.spectra.f_obs - self.sl_out.spectra.f_syn)
        self.axleft_low.plot(self.sl_out.spectra.l_obs, res, color='black')  # All
        aux_ = np.ma.masked_where(self.sl_out.spectra.f_wei < 0, res)
        self.axleft_low.plot(self.sl_out.spectra.l_obs, aux_, color='blue')  # If wei > 0
        aux_ = np.ma.masked_where(self.sl_out.spectra.f_wei != 0, res)
        self.axleft_low.plot(self.sl_out.spectra.l_obs, aux_, color='magenta')  # If wei == 0 -> Masked
        aux_ = np.ma.masked_where(self.sl_out.spectra.f_wei != -1, res)
        self.axleft_low.plot(self.sl_out.spectra.l_obs, aux_, 'x', color='red')  # If wei == -1 -> Clipped
        aux_ = np.ma.masked_where(self.sl_out.spectra.f_wei != -2, res)
        self.axleft_low.plot(self.sl_out.spectra.l_obs, aux_, 'o', color='green')  # If wei == -2 -> Flagged
        self.axleft_low.set_ylim(lim_res_low, lim_res_upp)

        ## RIGHT
        log_age = np.log10(self.sl_out.population.popage_base)
        Z = self.sl_out.population.popZ_base
        xj = self.sl_out.population.popx
        mu_cor = self.sl_out.population.popmu_cor

        ### Define bar spacing
        d = log_age[1:] - log_age[:-1]
        Zwidth = np.min(d[d > 0])

        ### RIGHT - UPPER
        aux_sum = np.zeros(np.shape(Z[Z == Z[0]]))
        i_color = 0
        for i_Z in np.unique(Z):
            v_ = (Z == i_Z)
            self.axright_upp.bar(log_age[v_], xj[v_], width=Zwidth, align='center', color=Zcolor[i_color],
                                 bottom=aux_sum, label=('%3.4f' % i_Z))
            aux_sum = aux_sum + xj[v_]
            i_color = i_color + 1

        ## RIGHT - LOWER
        aux_sum = aux_sum * 0.
        i_color = 0
        for i_Z in np.unique(Z):
            v_ = (Z == i_Z)
            aux_sum = aux_sum + mu_cor[v_]
            i_color = i_color + 1
        self.axright_low.bar(log_age[v_], aux_sum, width=Zwidth, align='center', color='white', label=('%3.4f' % i_Z),
                             log=True)
        self.axright_low.set_xlim(np.min(log_age) * .99, np.max(log_age) * 1.01)

        # Remove last lower ytick
        self.axleft_low.set_yticks(self.axleft_low.get_yticks()[:-1])
        self.axright_low.set_yticks(self.axright_low.get_yticks()[:-1])

        # Axis Labels
        self.axright_low.set_xlabel('log age [yr]')
        self.axright_low.set_ylabel('$\mu_j$ [%]')
        self.axright_upp.set_ylabel('$x_j$ [%]')

        self.axleft_low.set_xlabel('$\lambda [\AA]$')
        self.axleft_low.set_ylabel('Residual spectrum')
        self.axleft_upp.set_ylabel('$F_\lambda [normalized]$')

        # Some fit labels...
        self.fig.text(.6, .88, '$\chi^2 =\ $' + ('%3.2f' % self.sl_out.keywords['chi2']))
        self.fig.text(.6, .84, 'adev = ' + ('%3.2f' % self.sl_out.keywords['adev']))
        self.fig.text(.6, .80, '$S/N =\ $' + ('%3.2f' % self.sl_out.keywords['SN_normwin']))
        self.fig.text(.6, .76, '$A_V =\ $' + ('%3.2f' % self.sl_out.keywords['A_V']))
        self.fig.text(.6, .72,
                      '$\sigma_\star =\ ' + ('%3.2f' % self.sl_out.keywords['v_d']) + '\ $km/s\t$v_\star =\ ' + (
                                  '%3.2f' % self.sl_out.keywords['v_0']) + '\ $km/s')