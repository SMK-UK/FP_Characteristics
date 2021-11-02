"""
Sean Keenan, Heriot-Watt University, Edinburgh"
Mazzera Group Summer Project
Narrow Linewidth Laser Cavity Design
"""

# Import modules
from numpy import array, empty, linspace
from math import exp, pi, sqrt, cos, log
from numpy.polynomial.polynomial import polyfit, polyval
import matplotlib.pyplot as mp
from matplotlib import ticker

# set flag for which plots to generate, 0 for all, 1 for FSR and Linewidth, 2 for FP Spectrum
flag1 = 0
# set flag to save plots and path (0 no, 1 yes)
flag2 = 1

path = "C:\\Users\\sk88\\Documents\\Python\\FP_Cavity\\Images_2\\"
# speed of light in vacuum (m/s)
c = 2.9972e8
# reflectivity of mirror 1
r_1 = 0.9985
# reflectivity of mirror 2
r_2 = 0.9985
# square of r
R = sqrt(r_1 * r_2)
# absorption co-efficient
alpha = 0
# radius of curvature (m)
C = 200e-3
# mount length
l = 1e-3
# material expansion (invar and mount) 
dl = 15.42e-6
# cavity length (m)
L_0 = 250e-3
# expanded cavity length
L_1 = L_0 #L_0 + dl + 2*l
# resolution for length change and wavelength
res = [2e-6, 606*12, 1550*12]
# generate array of lengths around L_0
L = array([L_1 - res[0], L_1, L_1 + res[0]])
# desired lambda (m)
wave = array([606e-9, 1550e-9])
# number of data points for model
num_steps = 10000
# array of wavelengths for model
wave_array = array([linspace(start=wave[0] - 0.0006e-9 , stop=wave[0] + 0.0006e-9, num=num_steps),\
             linspace(start=wave[1] - 0.0015e-9, stop=wave[1] + 0.0015e-9, num=num_steps)])
# subsequent frequency (Hz)
nu_0 = c / wave
# decay time constant
tau_c = -(2 * L) / (c * log(R ** 2 * (1 - alpha) ** 2))
# output linewidth (Hz)
delta_nu = (2 * pi * tau_c) ** -1
# free spectral range (Hz)
FSR = c / (2 * L)
# finesse from FSR
finesse_0 = FSR / delta_nu
# finesse from reflectivity
finesse = pi * sqrt(R) / (1 - R)
# intensity (arb. units)
I_0 = 100
# effective index of medium
n_eff = 1

# create empty array for delta phase delay term
delta = empty([len(wave_array), len(L), num_steps], dtype=float)

# generate data for delta array
for m in range(len(wave_array)):
        for n in range(len(L)):
                for o in range(len(wave_array[m])):
                        delta[m][n][o] = (4 * pi * L[n]) / wave_array[m][o]

# create empty arrays for model
denominator = empty([len(L), num_steps], dtype=float)
I_t = empty([len(wave_array), len(L), num_steps], dtype=float)

# generate data for model
for p in range(len(wave_array)):
        for q in range(len(L)):
            numerator = I_0 * (1 - R * exp(-2 * alpha * L[q])) ** 2
            for r in range(len(wave_array[p])):
                denominator = 1 + R ** 2 * exp(-4 * alpha * L[q]) - 2 * R * exp(-2 * alpha * L[q]) * cos(delta[p][q][r])
                I_t[p][q][r] = numerator / denominator

# generate plots

# set global tick size
mp.rcParams['xtick.labelsize'] = 8
mp.rcParams['ytick.labelsize'] = 8

if flag1 ==0 or flag1 ==1:

        # interpolate data
        coeff_1 = polyfit(x=L, y=delta_nu, deg=2)
        coeff_2 = polyfit(x=L, y=FSR, deg=2)
        fit_1 = polyval(L, coeff_1)
        fit_2 = polyval(L, coeff_2)

        # create subplots and handles to axis
        fig_0, axs_0 = mp.subplots(nrows=1, ncols=2, sharex='all', sharey='none', constrained_layout='true')

        axs_0[0].plot(L * 1e3, delta_nu / 1e3, 'xr', markersize='8')
        axs_0[0].plot(L * 1e3, fit_1 / 1e3, color='red', linestyle='dashed', label='fit')
        axs_0[0].set_title('$\\Delta\\nu$ v. Length', fontsize='10')
        axs_0[0].set(ylabel='Frequency (KHz)')

        # plot data to axis 2 and format
        axs_0[1].plot(L * 1e3, FSR / 1e6, 'xb', markersize='8')
        axs_0[1].plot(L * 1e3, fit_2 / 1e6, color='blue', linestyle='dashed', label='fit')
        axs_0[1].set_title('FSR v. Length', fontsize='10')
        axs_0[1].set(ylabel='FSR (MHz)')

        # set formatting for both axis
        for ax in axs_0.flat:
                # plots gridlines
                ax.grid()
                # display legend
                ax.legend(loc='upper right')
                # set axis labels
                ax.set(xlabel='Length (mm)')
                # change tick label format
                ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))
                ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.3f}"))

        if flag2 ==1:
                fig_0.savefig(fname=path + 'LwidthFSR_v_Length.pdf', format='pdf')

if flag1 ==0 or flag1 ==2:

        fig_1, axs_1 = mp.subplots(nrows=2, ncols=1, sharex='none', sharey='none', constrained_layout='true')

        for s, length in enumerate(L):
                axs_1[0].plot(wave_array[0]*1e9, I_t[0][s][:], label='length'+str('{:.5f}'.format(L[s]*1e3))+'mm')
                axs_1[1].plot(wave_array[1]*1e9, I_t[1][s][:], label='length'+str('{:.5f}'.format(L[s]*1e3))+'mm')

        # set formatting for both axis
        for ax in axs_1.flat:
                # plots gridlines
                ax.grid()
                # display legend
                ax.legend()
                # set axis labels
                ax.set(xlabel='Wavelength (nm)')
                # set axis labels
                ax.set(ylabel='Intensity (A.U.)')
                # change tick label format
                ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4f}"))
                ax.xaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.4f}"))

        if flag2 ==1:
                fig_1.savefig(fname=path + 'FSR.pdf', format='pdf')

mp.show()