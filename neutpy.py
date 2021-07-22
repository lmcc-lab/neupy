import pandas as pd
import numpy as np
import math
import time
import constants
import preferences
import matplotlib.pyplot as plt
import plotResults
from dependencies import conv_str_to_list as ctl
from save_updated_data import contributers, maxContributions
import ContributingChartGenerator as ccg

if preferences.full == False:
    running_CDF, running_PDF, num_contrib, running_contributions, running_percent = maxContributions()
    ccg.gen_svg(running_contributions, running_percent)


def single_fission_event(time, update = False, plotResults = False):
    if update == True:
        update_data()
    CDFPDF_data = pd.read_csv('./Contributing_chains/CDF_PDF_full.csv')
    CDF = ctl(CDFPDF_data.iloc[0,1])
    PDF = ctl(CDFPDF_data.iloc[0,2])
    if plotResults == True:
        plt.plot(time, CDF, label = 'Full')
        if preferences.full == False:
            plt.plot(time, running_CDF, label = str('Select contributions'))
            plt.title('CDF with '+str(round(running_percent,2))+'% Chain contributions')
        else:
            plt.title('CDF')
        plt.xscale('log')
        plt.xlabel('time (s)')
        plt.grid()
        plt.legend()
        plt.ylabel('Cumulative neutrino emissions')
        plt.show()

        plt.plot(time, PDF, label = 'Full')
        if preferences.full == False:
            plt.plot(time, running_PDF, label = 'Select contributions')
            plt.title('PDF with '+str(round(running_percent,2))+'%'+' Chain contributions')
        else:
            plt.title('PDF')
        plt.xlabel('time (s)')
        plt.ylabel('Probability density of neutrino emissions')
        plt.xscale('log')
        plt.yscale('log')
        plt.legend()
        plt.grid()
        plt.show()
