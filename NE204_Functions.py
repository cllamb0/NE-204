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

import sklearn.cluster as sklclu
from sklearn.cluster import KMeans
from kneed import KneeLocator

import matplotlib.animation as animation
from matplotlib.animation import FFMpegWriter

from matplotlib.colors import Normalize
import matplotlib.cm as cm

try:
    import becquerel as bq
except:
    print('Becquerel not imported.\nNot necessary at all for functions in this file.')

plt.rcParams['figure.dpi'] = 300

source_colors = {'Co57': 'red', 'Co60': 'orange', 'Cs137': 'gold',
                 'Eu152': 'green', 'Th228': 'blue', 'Background': 'purple'}

try:
    os.mkdir('Plots')
    print('Created Plots directory')
except:
    print('Plots directory already exists')
    pass


"""
    Functions and Variables created during Lab 1
"""
def import_data(datafile, adc_max=16383):
    # Imports data into a numpy array based on whatever raw data is present in the h5 file
    f = h5py.File(datafile)

    data = np.array(f['raw_data'])
    f.close()

    sat_i = []
    for i in tqdm(range(len(data)), desc='Checking for Saturation', leave=False):
        if np.max(data[i]) == adc_max:
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

def gauss_fit(x, y, var=False):
    # Fits gaussian and returns fit parameters
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    perr = np.sqrt(np.diag(pcov))
    if var:
        return popt, perr
    return popt

def calibrate_energy(ch, a1, a2, a3, a4, a5, b):
    # Energy calibrate channel number given calibration constants
    return a1*ch + a2*ch**2 + a3*ch**3 + a4*ch**4 + a5*ch**5 + b

# Importing in pre-determined calibration constants
calib_consts = np.load('../Lab-1/calibration_values_final.npy')

def calibrate_pulses(data, peak=500, gap=2500, tau=15000, max_energy=3000, return_inds=False,
                    calibrated=True, max_trapezoid_height=1000, save_file=None, svw=51, use_CFD=False, shift=500,
                    samp_size=500, new_trapezoid=True, use_calib=None, multipulse_rejection=False, verbose=False):
    # Takes in raw_data array and creates the trapezoid, calculates the height
    # and then calibrates the energy corresponding to that trapezoid height
    if not calibrated:
        max_energy = max_trapezoid_height
    energies, rinds = [], []
    if multipulse_rejection:
        pileup_inds = possible_pileup(data, peak, gap, tau)
    else:
        pileup_inds = []
    for i, pulse in enumerate(tqdm(data, desc="Creating spectra", leave=False)):
        fs = savgol_filter(pulse, svw, 0)
        try:
            if i in pileup_inds:
                continue

            if new_trapezoid:
                if use_calib is None:
                    use_calib = calib_consts_new
                if i == 0:
                    if verbose:
                        print('Using new trapezoid filtering technique')
                trap = trapezoid_filter(fs, peak, gap, tau, samp_size=samp_size)
            else:
                if use_calib is None:
                    use_calib = calib_consts
                if i == 0:
                    if verbose:
                        print('Using old trapezoid filtering technique')

                if use_CFD:
                    trap = s(fs, CFD(fs, 0.1, shift=shift, samp_size=samp_size), tau, peak, gap)
                else:
                    trap = s(fs, determine_rise(fs), tau, peak, gap)
            #trap_height_fit = trapezoid_height(trap, peak, gap)
            if calibrated:
                pulse_height = calibrate_energy(np.max(trap), *use_calib)
            else:
                pulse_height = np.max(trap)
        except:
            if verbose:
                print('Index {} failed to fit or create trapezoid'.format(i))
            continue

        if pulse_height < max_energy:
            energies.append(pulse_height)

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

def quick_resolution(energies, e_min, e_max, adc_bit=11):
    # Goes through a quick gaussian fitting based on a calibrated energies array and min/max energies
    plt.figure(dpi=300)
    counts, bins_out, patches = plt.hist(energies, bins=2**adc_bit, alpha=0.5)

    win_min, win_max = e_min, e_max

    start = np.argmin(np.abs(bins_out - win_min))
    end = np.argmin(np.abs(bins_out - win_max))

    gauss_result, gauss_unc = gauss_fit(bins_out[start:end], counts[start:end], var=True)

    print('Peak energy: {} +- {}'.format(gauss_result[2], gauss_unc[2]))
    print('Peak sigma: {} +- {}'.format(gauss_result[3], gauss_unc[3]))
    print("Calculated resolution: {}%".format(round(np.abs(((gauss_result[3]*2.35482)/gauss_result[2])*100),3)))

    plt.plot(gauss(np.arange(3000), *gauss_result), alpha=0.5)
    plt.xlim(win_min, win_max)
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts')
    plt.semilogy()
    plt.show()

    return gauss_result, gauss_unc
"""
    Functions and Variables created during Lab 2
"""

def import_lab1_energies(source_name, indexes=False, new=True):
    if new:
        filepath_enes = '../Lab-1/Data/Trapezoid_Heights/New/'
    else:
        filepath_enes = '../Lab-1/Data/Trapezoid_Heights/Old/'
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

def CFD(signal, frac, shift=500, pf=False, samp_size=500, save_name=None):
    signal = reset_zero(signal, samp_size=samp_size)
    signal = signal[:np.argmax(signal)+shift]# For cropping before decay

    mod_signal = (-frac*signal)[shift:]
    CFD_result = signal[:len(signal)-shift]+mod_signal

    if pf:
        plt.plot(4e-3*np.arange(len(signal)), signal, label='Original Signal')
        plt.plot(4e-3*np.arange(len(mod_signal)), mod_signal, label='Modified Signal')
        plt.plot(4e-3*np.arange(len(CFD_result)), CFD_result, label='CFD Signal')
        plt.axvline(4e-3*np.argwhere(np.where(CFD_result>0, 0, CFD_result)!=0)[-1][0], label='CFD {} Result'.format(frac), color='red')
        plt.xlabel('Time [us]')
        plt.ylabel('Pulse Height [arb]')
        plt.legend(loc='lower right')
        if save_name is not None:
            plt.savefig('Plots/{}.png'.format(save_name), dpi=300, facecolor='white', bbox_inches='tight')
        #return mod_signal, CFD_result, np.argwhere(np.where(CFD_result>0, 0, CFD_result)!=0)[-1][0]
    else:
        return np.argwhere(np.where(CFD_result>0, 0, CFD_result)!=0)[-1][0]

def tangent_line(x, y, point, width=10):
    ind = x.index(point)
    mb = (y[ind+1]-y[ind])/(x[ind+1]-x[ind])
    ma = (y[ind]-y[ind-1])/(x[ind]-x[ind-1])
    m = (mb+ma)/2

    xt = np.linspace(point-width/2, point+width/2, 1000)
    yt = m*(xt - x[ind]) + y[ind]

    return xt, yt

def cluster_data(data, n_clust, buff=50, quality=True):
    kmeans = sklclu.KMeans(init='random', n_clusters=n_clust, n_init=10, max_iter=2000, random_state=69) #int(time.time()))
    kmeans.fit(data)

    clusts = kmeans.labels_

    rise_times = []
    for pp in range(len(data)):
        signal = savgol_filter(data[pp], 31, 0)
        rise_times.append([CFD(signal, 0.1, samp_size=50), CFD(signal, 0.9, samp_size=50)])

    # Attempting FOM creation
    if quality:
        w_inertia = []
        for i in range(len(clusts)):
            weight = np.zeros(len(data[i]))
            weight[rise_times[i][0]:rise_times[i][1]+buff] = np.ones((rise_times[i][1]+buff)-(rise_times[i][0]))
            weight = weight/sum(weight)
            #val = weight*np.abs((data[i]-kmeans.cluster_centers_[clusts[i]])/kmeans.cluster_centers_[clusts[i]])
            #sse = np.cumsum(np.square(data[i]-kmeans.cluster_centers_[kmeans.labels_[i]])/kmeans.cluster_centers_[kmeans.labels_[i]])[-1]
            val = np.cumsum(weight*np.square(data[i]-kmeans.cluster_centers_[clusts[i]]))[-1]
            w_inertia.append(val)

        return np.sum(np.array(w_inertia)), kmeans.inertia_
    else:
        return kmeans

### The return of the trapezoids and my dissatifaction for the original trapezoid filtering that I made

def MWD(pulse, M, tau, samp_size=500):
    # Width of window M must be greater than the rise time of the pulse
    signal = delay_signal(reset_zero(pulse, samp_size=samp_size), M)
    sum_sig = np.cumsum(signal)
    return signal[M:]-signal[:-M]+(1/tau)*(sum_sig[M:]-sum_sig[:-M])

def MA(pulse, L, samp_size=500):
    signal = delay_signal(reset_zero(pulse, samp_size=samp_size), L)
    sum_sig = np.cumsum(signal)
    return (1/L)*(sum_sig[L:]-sum_sig[:-L])

def trapezoid_filter(pulse, peaking_time=500, gap_time=2500, tau=15000, samp_size=500):
    return MA(MWD(pulse, peaking_time+gap_time, tau, samp_size=samp_size), peaking_time, samp_size=samp_size)

def plot_spectra(energies, names=None, calibrated=True, adc_bit=10):
    """
    If energies is an array of pulse heights or energies, this will plot a single spectra
    If energies is a list of array of phs or energies then this will plot multiple spectra given "names"
    """
    if type(energies[0]) == np.ndarray:
        bins_ = np.linspace(0, np.max(energies[0]), 2**adc_bit+1)
        if names is None:
            names = ['Spectra {}'.format(i) for i in np.arange(len(energies))+1]
        for e in range(len(energies)):
            plt.hist(energies[e], bins=bins_, label=names[e], alpha=0.5)
    else:
        bins_ = np.linspace(0, np.max(energies), 2**adc_bit+1)
        if names is None:
            names = 'Spectra'
        plt.hist(energies, bins=bins_, label=names)

    plt.semilogy()
    plt.legend()
    if calibrated:
        plt.xlabel('Energy [keV]')
    else:
        plt.xlabel('Trapezoid Height')
    plt.ylabel('Counts')

# Importing in pre-determined calibration constants from improved trapezoidal filtering
calib_consts_new = np.load('../Lab-2/calibration_values_new.npy')

def possible_pileup(pulses, peak, gap, tau, fudge=1.05, fudge_max=3000, binary=False):
    # Determine indices where there is possible pileup or multiple gamma rays
    pileup_inds = []
    if binary:
        theo_area = np.max(reset_zero(pulses))*(peak+gap)
        trap = trapezoid_filter(pulses, peak, gap, tau)

        if (np.sum(trap) > fudge*theo_area) or (np.argmax(pulses) > fudge_max):
            return True
        return False

    for p, pulse in enumerate(tqdm(pulses, leave=False, desc='Detecting Pileup')):
        theo_area = np.max(reset_zero(pulse))*(peak+gap)
        trap = trapezoid_filter(pulse, peak, gap, tau)

        if (np.sum(trap) > fudge*theo_area) or (np.argmax(pulse) > fudge_max):
            pileup_inds.append(p)

    return pileup_inds

def better_calibrate_pulses(data, kmeans, optimums, tau=15000, return_inds=False, svw=51,
                            samp_size=500, max_energy=3000, use_calib=calib_consts_new, save_file=None,
                            shifting=True, old=False, ignore_clusters=[]):
    # Optimums should be a nested list containing the optimum shaping parameters and energy shifts
    # for each of the n clusters within kmeans where each element: [[peak, gap], shift]

    mod_data, mod_inds = [], []
    for p in tqdm(range(len(data)), desc='Modifying data for kmeans identification', leave=False):
        try:
            start = determine_rise(savgol_filter(data[p][:2200], 51, 0))
            sig = reset_zero(data[p][start:start+1100], samp_size=50)
            mod_data.append(sig/max(sig))
            mod_inds.append(p)
        except:
            print('Index {} failed for some reason?'.format(p))
            pass
    mod_inds = np.array(mod_inds)
    mod_data = np.array(mod_data)

    mod_clusters = kmeans.predict(mod_data)

    end_inds, pileup_inds, ignored_inds = [], [], []
    energies = []
    for mp, pulse_ind in enumerate(tqdm(mod_inds, desc='Creating better spectra', leave=False)):
        cl = mod_clusters[mp]
        if cl in ignore_clusters:
            ignored_inds.append(pulse_ind)
            continue
        # Checking for pileup
        if old:
            peak, gap = 500, 2500
        else:
            peak, gap = optimums[cl][0][0], optimums[cl][0][1]

        if possible_pileup(data[pulse_ind], peak, gap, tau, binary=True):
            pileup_inds.append(pulse_ind)
            continue

        trap = trapezoid_filter(data[pulse_ind], peak, gap, tau, samp_size=samp_size)

        if shifting:
            energy = calibrate_energy(np.max(trap), *use_calib) + optimums[cl][1]
        else:
            energy = calibrate_energy(np.max(trap), *use_calib)

        if energy <= max_energy:
            energies.append(energy)
            if return_inds:
                end_inds.append(pulse_ind)

    energies = np.array(energies)

    if save_file is not None:
        try:
            os.mkdir('Data/Better_Spectra/')
        except:
            pass
        print('Saving better energies')
        np.save('Data/Better_Spectra/'+save_file+'.npy', energies)
        if return_inds:
            np.save('Data/Better_Spectra/'+save_file+'-indexes.npy', end_inds)

    if return_inds:
        end_inds = np.array(end_inds)
        pileup_inds = np.array(pileup_inds)
        ignored_inds = np.array(ignored_inds)
        return energies, end_inds, pileup_inds, ignored_inds
    return energies

def combine_data_files(input_dir, input_file_hint, output_directory):
    """
        input_file_hint: some text that shows up in all the files to combine
    """
    uncombined_data = [input_dir+f for f in os.listdir(input_dir) if f != '.DS_Store' and input_file_hint in f]

    file_name = '-'.join(uncombined_data[0].split('/')[-1].split('_')[0].split('-')[:-1])

    datas = []
    for i, file in enumerate(uncombined_data):
        f = h5py.File(file, 'r')
        data_np = np.array(f['raw_data'])
        time_event, unique_times = np.unique(np.array(f['event_data']['timestamp']), return_index=True)

        data_np_crop = data_np[unique_times]

        datas.append(data_np_crop)
    merge_data = np.concatenate(datas)

    fmerged = h5py.File('{}/{}-{}.h5'.format(output_directory, file_name, i+1), 'w')
    fmerged.create_dataset('raw_data', data=merge_data)

    fmerged.close()

    print('Combined {} files of data'.format(i+1))


def optimize_trapezoid_parameters(data, emin, emax, svw=31):
    # Optimizes peaking and gap time for a set of pulses given an energy range of interest
    # Used given purely photopeak events within a spectra
    rise_times = []
    for pp in tqdm(range(len(data)), desc='Determining rise times in data', leave=False):
        signal = savgol_filter(data[pp], svw, 0)
        rise_times.append(CFD(signal, 0.9, samp_size=50)-CFD(signal, 0.1, samp_size=50))
    rise_times = np.array(rise_times)

    peak_range = range(200, 610, 10)
    gap_range = range(500, 2510, 10)

    bins_ = np.linspace(emin, emax, (emax-emin)+1)
    means, stds = [], []
    parameters = []
    for p in tqdm(peak_range, leave=False):
        means_p, stds_p = [], []
        for g in tqdm(gap_range, leave=False):
            energies_temp = calibrate_pulses(data, p, g, 15000)
            counts, b, _ = plt.hist(energies_temp, bins=bins_)
            plt.close()
            gf = gauss_fit(bins_[:-1], counts)

            means_p.append(gf[2])
            stds_p.append(gf[3])
            parameters.append([p, g])
        means.append(means_p)
        stds.append(stds_p)
    means = np.array(means)
    stds = np.array(stds)
    return means, stds, parameters
