# Import libraries

import os
import traitlets
from ipywidgets import *
from IPython.display import display, Javascript, clear_output

import numpy as np
import pandas as pd
import matplotlib
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline                    # Safe smooth interpolation
from scipy.signal import find_peaks
from scipy.signal import medfilt
from scipy import signal
from pybaselines import Baseline, utils
from scipy.integrate import simpson, trapz

from tkinter import Tk, filedialog

import peakutils
from peakutils.plot import plot as pplot

from tqdm import tqdm

## Defining all functions

# DO NOT CHANGE anything - except the locaiton of your lamp data file ( it is in the i_corr function)

class SelectFilesButton(widgets.Button):
    """A file widget that leverages tkinter.filedialog."""

    def __init__(self, *args, **kwargs):
        """Initialize the SelectFilesButton class."""
        super(SelectFilesButton, self).__init__(*args, **kwargs)
        # Add the selected_files trait
        self.add_traits(files=traitlets.traitlets.List())
        # Create the button.
        self.description = "Select File"
        self.icon = "square-o"
        self.style.button_color = "orange"
        # Set on click behavior.
        self.on_click(self.select_files)

    @staticmethod
    def select_files(b):
        """Generate instance of tkinter.filedialog.
        Parameters
        ----------
        b : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        # List of selected fileswill be set to b.value
        b.files = filedialog.askopenfilename(multiple=True)
        f = b.files
        '''
        for i in range(0,len(f)):
            print(i,f[i],end="\n")   # load preview cut process plot
        '''
        b.description = "Files Selected"
        b.icon = "check-square-o"
        b.style.button_color = "lightgreen"
 
class DataReader():
    
    def __init__(self, file_name):
        self.file_name = file_name
        
    def read_file(self):
        
        fname = os.path.splitext(self.file_name)[0]
        ext = os.path.splitext(self.file_name)[1]
        
        ext_dict = {".csv" : (",",0), ".txt" : ("\t",1), ".dat" : (None,0)}
 
        data = pd.read_csv(self.file_name, skiprows=ext_dict[ext][1], sep=ext_dict[ext][0],engine='python',header=None)
        x = data[0].to_numpy()
        y = data[1].to_numpy()
        
        if ext == ".dat":
            col_len = len(data.columns)
            
            y = np.zeros(shape=(len(x),col_len-1))

            for k in range(0,col_len-1):
                y[:,k] = data.iloc[:,k+1].values
        
        return x, y

    def read_xlsx(self):
        return pd.read_excel(self.file_name)

    def read_json(self):
        return pd.read_json(self.file_name)

def data_cleaning(data,erp):
    
    # Code to remove cosmic rays
    #erp=100 # Set percentage error limit
    err=erp/100

    Ndatapts = len(data)
    ydata = data

    for i in range(0,Ndatapts): # This loop removes hot regions one pixel wide
        if i>0 and i<Ndatapts-1:
            if ydata[i]>(1+err)*ydata[i+1] and ydata[i]>(1+err)*ydata[i-1]:
                #print('Entered the loop 1')
                ydata[i]=(ydata[i-1]+ydata[i+1])/2

    for i in range(0,Ndatapts): # This loop removes hot regions up to 3 pixels wide
        if i>1 and i<Ndatapts-2 and ydata[i]>(1+err)*ydata[i+2] and ydata[i]>(1+err)*ydata[i-2]:
            #print('Entered the loop 3')
            ydata[i], ydata[i-1], ydata[i+1] = (ydata[i-2]+ydata[i+2])/2, (ydata[i]+ydata[i-2])/2, (ydata[i+2]+ydata[i])/2

    for i in range(0,Ndatapts): # This loop removes hot regions up to 5 pixels wide
        if i>4 and i<Ndatapts-5 and ydata[i]>(1+err)*ydata[i+5] and ydata[i]>(1+err)*ydata[i-5]:
            #print('Entered the loop 3')
            ydata[i], ydata[i-1], ydata[i-2], ydata[i-3], ydata[i-4], ydata[i+1], ydata[i+2], ydata[i+3], ydata[i+4] = (ydata[i-5]+ydata[i+5])/2, (ydata[i]+ydata[i-2])/2, (ydata[i-1]+ydata[i-3])/2, (ydata[i-2]+ydata[i-4])/2, (ydata[i-3]+ydata[i-5])/2, (ydata[i+2]+ydata[i])/2, (ydata[i+3]+ydata[i+1])/2, (ydata[i+4]+ydata[i+2])/2, (ydata[i+5]+ydata[i+3])/2

    for i in range(0,Ndatapts): # This loop removes dead regions one pixel wide
        if i>0 and i<Ndatapts-1 and ydata[i]<(1-err)*ydata[i+1] and ydata[i]<(1-err)*ydata[i-1]:
            #print('Entered the loop 4')
            ydata[i]=(ydata[i-1]+ydata[i+1])/2

    for i in range(0,Ndatapts): # This loop removes dead regions up to 3 pixels wide
        if i>1 and i<Ndatapts-2 and ydata[i]<(1-err)*ydata[i+2] and ydata[i]<(1-err)*ydata[i-2]:
            #print('Entered the loop 5')
            ydata[i], ydata[i-1], ydata[i+1] =(ydata[i-2]+ydata[i+2])/2, (ydata[i]+ydata[i-2])/2, (ydata[i+2]+ydata[i])/2
            
    return ydata

def subtract_bsl(x,y,hf):
    
    baseline_fitter = Baseline(x, check_finite=False)
    base = baseline_fitter.mor(y, half_window=hf)[0]
    
    return base

def filter_savgol(data,k_size,p_order):
    
    data[np.isnan(data)] = 0
    
    sav_filt_data=signal.savgol_filter(data,k_size,p_order)
    
    return sav_filt_data

def filter_median(data,k_size):
    
    filt_data = medfilt(data,kernel_size=k_size)
    
    return filt_data

def i_corr(flamp,f,i):
    start = 0
    end = -1

    head_i, tail_i = os.path.split(f[i])
    
    ## Important note to change the location of this .txt file below to where it is in your computer.
    calibstd = np.loadtxt(r"G:\Shared drives\Pauzauskie Team Drive\Users\CG\Scripts\030410638_HL-2000-CAL_2014-01-15-14-09_VISEXT1EXT2_FIB.txt")
    xcalib = calibstd[:,0]
    ycalib = calibstd[:,1]

    x,y = DataReader(file_name=f[i]).read_file()

    HglampFunc = CubicSpline(xcalib,ycalib)
    hglampI = HglampFunc(x) # Create interpolation of true lamp spectrum
    
    hglampdata_x, hglampdata_y = DataReader(file_name=flamp[0]).read_file() # Split true lamp spectra into x and y

    ICF = hglampI/(hglampdata_y) # Creates ratio of true lamp spectra to real lamp data, ICF = Intensity Correction Factor

    ynew = (y)*ICF # multiplies real data by intensity correction factor

    ynew = np.nan_to_num(ynew,nan=0,posinf=0,neginf=0)
    datamatrix = np.column_stack((x,ynew)) # Compiles corrected data into a new matrix
    savename = f[i][:-4]+"_calib.txt" # Create filename for new data
    
    head_c, tail_c = os.path.split(savename)
    
    np.savetxt(savename, datamatrix) # Save new data

    return x,y,ynew

def slicing_func(x,y,yi,start_x,end_x):
    start = np.argmin(abs(x-start_x))
    end = np.argmin(abs(x-end_x))
    #print(start,end)
    x_cut = x[start:end]
    y_cut = y[start:end]
    yi_cut = yi[start:end]
    
    return x_cut, y_cut, yi_cut, start, end

def _3Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2, amp3,cen3,wid3):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
            (amp2*wid2**2/((x-cen2)**2+wid2**2)) +\
                (amp3*wid3**2/((x-cen3)**2+wid3**2))

def norm_func(x, y, do_normalize):
    
    if do_normalize == 'relative':
        norm_intens_data = (y - np.min(y)) / (np.max(y) - np.min(y))
    elif do_normalize == 'scale':
        norm_intens_data = y / np.max(y)
    elif do_normalize == 'area':
        norm_intens_data = y/np.trapz(y,x)
    else:
        norm_intens_data = y
    
    return norm_intens_data

def mean_f_wvl(x,y,meanf_method):
    
    lambda_f = None
    
    if meanf_method == 'area':
        lambda_f = np.sum(x*y**2) /np.sum(y**2)
    elif meanf_method == 'sum':
        lambda_f = np.sum(x*y)/np.sum(x)
    elif meanf_method == 'integrate':
        lambda_f = simpson(x*y,x)/simpson(y,x)
    elif meanf_method == 'energy':
        x = 1240/x
        y = (1240/(x**2))*y
        e_f = simpson(x*y,x)/simpson(y,x)
        lambda_f = 1240/e_f
        
    return lambda_f
                
def _1Lorentzian(x, amp1, cen1, wid1):
    return (amp1*wid1**2/((x-cen1)**2+wid1**2))

def ram2nm(x,laser_source):
    return 10**7/(10**7/laser_source - x)

def nm2ram(x,laser_source):
    return 10**7*(1/laser_source - 1/x)

def plot_meanf_data(laser_Pow, meanf_df, ax2, plot_title):
    for laser_power in laser_Pow:
        #power_str = str(laser_power)
        df2 = meanf_df.iloc[list(meanf_df["Power"]==laser_power)]
        meanf_val = df2['Mean Wvl'].mean()
        sme_val = df2['Mean Wvl'].std()/df2['Mean Wvl'].count()
        ax2.errorbar(laser_power,meanf_val,yerr=sme_val, marker='s',color='k',markersize=10,capsize=6,barsabove=True,ecolor='k',elinewidth=2,capthick=1)
        
    ax2.set_title(plot_title)
    ax2.set_xlabel('1020 nm laser power')
    ax2.set_ylabel('Mean f wavelength')
    #ax2.set_xlim()
    #ax2.set_ylim()
    
    #ax2.set_yticks([])
    #ax2.tick_params(axis='y',which='both',left='off',labelleft='off')

    #ax2.legend(bbox_to_anchor=(1, 1),title='')
    
    return ax2


def i_corr_cleaning(iFiles,do_baseline_subtraction,do_median_filtering,do_data_cleaning,do_savgol_filtering,do_conversion,do_normalize,do_peak_finding,p_order,k_size,erp,k_size_savgol,p_order_savgol):
    idata = iFiles.files

    for i in tqdm(range(len(idata))):
        
        idata_name = idata[i]
        print(idata_name)

        head_i, tail_i = os.path.split(idata[i])

        ixdata, iydata = DataReader(file_name=idata[i]).read_file()

        if do_median_filtering == 'y':
            filt_idata = filter_median(iydata, k_size)
        else:
            filt_idata = iydata

        if do_baseline_subtraction == 'y':
            base = subtract_bsl(filt_idata, p_order)
        else:
            base = 0

        bsl_subt_idata = filt_idata - base

        if do_data_cleaning == 'y':
            cleaned_idata = data_cleaning(bsl_subt_idata, erp)
        else:
            cleaned_idata = bsl_subt_idata

        if do_savgol_filtering == 'y':
            savgol_filt_idata = filter_savgol(cleaned_idata, k_size_savgol,p_order_savgol)
        else:
            savgol_filt_idata = cleaned_idata

        ixval, iyval = ixdata, savgol_filt_idata

        iyval_new = np.nan_to_num(iyval,nan=0,posinf=0,neginf=0)
        i_datamatrix = np.column_stack((ixval,iyval_new)) # Compiles corrected data into a new matrix
        i_save = idata_name[:-4]+"_cleaned.txt" # Create filename for new data
        print(i_save)
        head_is, tail_is = os.path.split(i_save)
        np.savetxt(i_save, i_datamatrix) # Save new data
        
        #ax1.plot(ixval, iyval_new,label=tail_is) # Code line to plot the curves
    # ax1.set_title("Intensity Correction files")
    # ax1.set_xlabel(x_label)
    #ax1.set_ylabel('Intensity (counts)')

    # ax1.set_xlim(3000,4000)
    # ax1.set_ylim(1.5e-5,3.5e-5)

    # ax1.legend(bbox_to_anchor=(0.8, 0.8), loc='center', title="File Name")

    # ax1.grid(color = 'green', which='both', linestyle = '-.', linewidth = 0.5, alpha=0.4)
        
    # plt.show()
    # plt.savefig()


if __name__ == "__main__":
    pass