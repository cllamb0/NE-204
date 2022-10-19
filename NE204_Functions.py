"""
    Bunch of functions created by Chris Lamb
    for the various labs in NE 204
"""

# Package Importing
import numpy as np
import matplotlib.pyplot as plt
import h5py
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
from tqdm.notebook import tqdm
import math
import warnings
warnings.filterwarnings("ignore")
import time
import os
from datetime import date

try:
    os.mkdir('Plots')
    print('Created Plots directory')
except:
    print('Plots directory already exists')
    pass


"""
    Functions and Variables created during Lab 1
"""
def import_data(datafile):
    # Imports data into a numpy array based on whatever raw data is present in the h5 file
    f = h5py.File(datafile)

    data = np.array(f['raw_data'])
    f.close()

    sat_i = []
    for i in range(len(data)):
        if data[i][1300] == 16383:
            sat_i.append(i)

    return data, sat_i

def determine_rise(signal, sigma=8, window=20, offset=100):
    # Input filtered signal, returns index of start of rise on pulse
    noise_samp = signal[:500]
    mean, std = np.mean(noise_samp), np.std(noise_samp)

    grad = np.gradient(np.where(signal > mean+sigma*std, signal, 0))

    grad_pos, grad_neg = np.argwhere(grad>2), np.argwhere(grad<-2)

    rise_start = 0
    for gp in grad_pos:
        close = False
        for gn in grad_neg:
            if gn-gp < window and gn-gp > 0:
                close = True
        if not close:
            rise_start = gp
            break

    return int(rise_start-offset)

def delay_signal(signal, delay=2000, sample_size=500, seed=69):
    # Adds an artifcial delay on pulse to act like 3316 delay functionality
    np.random.seed(seed)
    noise_samp = signal[:sample_size]
    mean, std = np.mean(noise_samp), np.std(noise_samp)
    noise = np.random.normal(mean, 0.9*std, delay)
    return np.hstack([noise, signal])

def dkl(signal, i, k, l, w):
    """
    Calculates dkl value given a signal,
        i = index of start of rise
        k = peaking time
        l = peaking time + gap time
        w = width of window to sample on pulse (trapezoid output width)
    """
    signal = delay_signal(signal, delay=w)
    vj = signal[i+w:i+w+w]
    vjk = signal[i-k+w:i+w+w-k]
    vjl = signal[i-l+w:i+w+w-l]
    vjkl = signal[i-k-l+w:i+w+w-k-l]
    return vj - vjk - vjl + vjkl

def s(signal, start_rise, tau, peaking_time, gap_time, w=0):
    """
    Input exponential signal, returns "trapezoid"
        start_rise = index of start of rise
        tau = decay constant of exponential (RC circuit related)
        peaking_time = "width" of rise of trapezoid
        gap_time = "width" of flat top of trapezoid
        w = width of trapezoid output width (default value calculated)
    """
    if w == 0:
        w = int(round(3*peaking_time+1.5*gap_time, 0))
    dkl_s = dkl(signal, start_rise, peaking_time, peaking_time+gap_time, w)

# New double big sigma method of calulating s below
#     dkl_sum, dkl_tau = np.cumsum(dkl_s), dkl_s*tau
#     return np.cumsum(np.array([dkl_sum[i] + dkl_tau[i] for i in range(w)]))
# Old recursive method of calculating s below
    ss = []
    for j in range(w):
        if j == 0:
            ss.append(0)
        else:
            ss.append(ss[j-1]*(1+1/tau)+ dkl_s[j])
    return np.array(ss)

def trapezoid(x, top_left, top_right, top, ps, pe, e):
    # Trapezoid signal shape
    tr = np.zeros(len(x))
    tr[:int(ps)] = 0
    tr[int(ps):int(top_left)] = ((top-0)/(top_left-ps))*(x[int(ps):int(top_left)]-ps)+(0)
    tr[int(top_left):int(top_right)] =  top
    tr[int(top_right):int(pe)] = ((e-top)/(pe-top_right))*(x[int(top_right):int(pe)]-pe) + e
    tr[int(pe):] = e
    return tr

def trapezoid_height(trapezoid_signal, peak=500, gap=500):
    # Fits trapezoid signal and then returns the average height of the flat top
    popt, pcov = curve_fit(trapezoid, np.arange(len(trapezoid_signal)), trapezoid_signal, method = 'lm',
                      p0=[peak, peak+gap, max(trapezoid_signal), 1/4*peak, len(trapezoid_signal)-200, 0])
    trap_height = np.mean(trapezoid_signal[math.floor(popt[0]):math.ceil(popt[1])])
    return trap_height

def gauss(x, H, A, x0, sigma):
    # Gaussian signal shape
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_fit(x, y):
    # Fits gaussian and returns fit parameters
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    perr = np.sqrt(np.diag(pcov))
    return popt

def calibrate_energy(ch, a1, a2, a3, a4, a5, b):
    # Energy calibrate channel number given calibration constants
    return a1*ch + a2*ch**2 + a3*ch**3 + a4*ch**4 + a5*ch**5 + b

# Importing in pre-determined calibration constants
calib_consts = np.load('../Lab-1/calibration_values_final.npy')

def calibrate_pulses(data, peak=500, gap=2500, tau=15000, max_energy=3000, return_inds=False,
                    calibrated=True, max_trapezoid_height=493154, save_file=None):
    # Takes in raw_data array and creates the trapezoid, calculates the height
    # and then calibrates the energy corresponding to that trapezoid height
    if not calibrated:
        max_energy = max_trapezoid_height
    energies, rinds = [], []
    for i, pulse in enumerate(tqdm(data, desc="Creating spectra")):
        fs = savgol_filter(pulse, 51, 0)
        try:
            trap = s(fs, determine_rise(fs), tau, peak, gap)
            #trap_height_fit = trapezoid_height(trap, peak, gap)
            if calibrated:
                pulse_energy = calibrate_energy(max(trap), *calib_consts)
            else:
                pulse_energy = max(trap)
        except:
            print('Index {} failed to fit or create trapezoid'.format(i))
            continue

        if pulse_energy < max_energy:
            energies.append(pulse_energy)

            if return_inds:
                rinds.append(i)

    energies = np.array(energies)

    if save_file is not None:
        try:
            os.mkdir('Data/Trapezoid_Heights/')
        except:
            pass
        print('Saving trapezoid heights')
        np.save('Data/Trapezoid_Heights/'+save_file+'.npy', energies)

    if return_inds:
        rinds = np.array(rinds)
        return energies, rinds
    return energies

def calc_activity(source):
    # Takes in info from "source_info" dictionary
    # [Production_Date, Measurement_Date, Half-Life(years), Produced_Activity]
    # returns activity at measurement date in Bq
    prod = date(*np.roll(np.array(list(map(int, source[0].split('/')))),1))
    meas = date(*np.roll(np.array(list(map(int, source[1].split('/')))),1))
    delt = ((meas-prod).days)/365.25

    return source[3]*np.exp(-(np.log(2)/source[2])*delt)*3.7e10

def quick_resolution(energies, e_min, e_max):
    # Goes through a quick gaussian fitting based on a calibrated energies array and min/max energies
    plt.figure(dpi=300)
    counts, bins_out, patches = plt.hist(energies, bins=2**11, alpha=0.5, label='Max')

    win_min, win_max = e_min, e_max

    start = np.argmin(np.abs(bins_out - win_min))
    end = np.argmin(np.abs(bins_out - win_max))

    gauss_result = gauss_fit(bins_out[start:end], counts[start:end])

    print('Peak energy: {}'.format(gauss_result[2]))
    print("Calculated resolution: {}%".format(round(np.abs(((gauss_result[3]*2.35482)/gauss_result[2])*100),3)))

    plt.plot(gauss(np.arange(3000), *gauss_result), alpha=0.5)
    plt.xlim(win_min, win_max)
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.semilogy()
    plt.show()


"""
    Functions and Variables created during Lab 2
"""

def import_lab1_energies(source_name, indexes=False):
    filepath_enes = '../Lab-1/Data/Trapezoid_Heights/'
    files_enes = ['Co57-Cal.npy', 'Co60-Cal.npy', 'Cs137-Cal.npy', 'Eu152-Cal.npy', 'Th228-Cal.npy', 'Background-Cal.npy']
    files_enes = [filepath_enes+f for f in files_enes]

    source_datasets = ['Co57', 'Co60', 'Cs137', 'Eu152', 'Th228', 'Background']
    if source_name in source_datasets:
        print('Loading in calibrated {} spectra'.format(source_name))
        energies = np.load(files_enes[source_datasets.index(source_name)])
        if indexes:
            energies_indexes = np.load(filepath_enes+'{}-indexes.npy'.format(source_name))
            return energies, energies_indexes
        return energies
    else:
        print('Enter valid lab 1 dataset from:\nCo57, Co60, Cs137, Eu152, Th228, Background')

def reset_zero(signal, samp_size=500):
    return signal - np.mean(signal[:samp_size])

def CFD(signal, frac, shift=500, pf=False, samp_size=500):
    signal = reset_zero(signal, samp_size=samp_size)
    mod_signal = (-frac*signal)[shift:]
    CFD_result = signal[:len(signal)-shift]+mod_signal

    if pf:
        return mod_signal, CFD_result, np.argwhere(np.where(CFD_result>0, 0, CFD_result)!=0)[-1][0]
    else:
        return np.argwhere(np.where(CFD_result>0, 0, CFD_result)!=0)[-1][0]
