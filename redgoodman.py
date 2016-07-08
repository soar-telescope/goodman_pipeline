#!/usr/bin/python3.4
#Simon Torres 2016-06-28
#pipeline for GOODMAN spectra reduction.



import sys
import os
import glob
#import pyfits as fits
import numpy as np
#import scipy.stats as stats
import re
#import subprocess
import time
import astropy.stats as asst
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')
#cython stuff
#import pyximport; pyximport.install()


#"global" variables
patt_over_trim  = 'to_'
patt_bias_corr  = 'bc_'
patt_flat_corr  = 'fc_'


class night:

    def __init__(self,date):
        self.date    = date
        self.all     = []
        self.bias    = []
        self.flats_wg= []
        self.flats_ng= []
        self.sci     = []
        self.lamp    = []
        self.master_bias = ''
        self.master_flat = ''
        self.last_prefix = ''
    
    #methods to store data
    def add_bias(self,inbias):
        self.bias.append(inbias)
        self.all.append(inbias)
        
    def add_flat_wg(self,inflat):
        self.flats_wg.append(inflat)
        self.all.append(inflat)
        
    def add_flat_ng(self,inflat):
        self.flats_ng.append(inflat)
        self.all.append(inflat)
        
    def add_sci(self,insci):
        self.sci.append(insci)
        self.all.append(insci)
    
    def add_lamp(self,inlamp):
        self.lamp.append(inlamp)
        self.all.append(inlamp)
    
    def set_master_bias(self,master):
        self.master_bias = master
     
    def set_master_flat(self,master):
        self.master_flat = master
        
    def set_last_prefix(prefix):
        self.last_prefix = prefix
        
    #methods to retrieve data 
    def get_bias(self):
        return self.bias
    
    def get_flats_wg(self):
        return self.flats_wg
    
    def get_flats_ng(self):
        return self.flats_ng

    def get_sci(self):
        return self.sci 
    
    def get_lamp(self):
        return self.lamp
    
    def get_all(self):
        return self.all
    
    
    
    
    

def get_file_list():
    lista = sorted(glob.glob("0*.fits")) #will need work
    try:
        header0 = fits.getheader(lista[0])
        this_night = night(header0["DATE"])
        for i in lista:
            #print("./" +i)
            header  = fits.getheader(i)
            obstype = header["OBSTYPE"]
            if obstype == "BIAS":
                this_night.add_bias(i)
            elif obstype == "FLAT":
                grating = header["GRATING"]
                if "NO GRATING" in grating:
                    this_night.add_flat_ng(i)
                elif grating == "KOSI_600":
                    this_night.add_flat_wg(i)
                else:
                    print("Unrecognized grating")
            elif obstype == "OBJECT":
                this_night.add_sci(i)
            elif obstype == "COMP":
                this_night.add_lamp(i)
        return this_night
    except IOError as err:
        if str(err) == "Empty or corrupt FITS file":
            print("Raised an IOError as ",err)
            os.system("rem.csh")
            sys.exit("Please run this script again")
        else:
            print(err)
            sys.exit("Please correct the errors and try again")
        
#GET DATA
def get_data_header(file_name):
  try:
    hdu0	= fits.open(file_name)
    scidata	= hdu0[0].data
    header	= hdu0[0].header
    scidata = scidata.byteswap().newbyteorder().astype('float64')
    return scidata,header
  except IOError as err:
    debug("I/O error (get_data_header): %s"%err)
  except TypeError as err:
    debug("TypeError (get_data_header): %s : %s"%(file_name,err))
    
    
def overscan_trim_corr(list_all):
    print_spacers("overscan and trim")
    length = len(list_all)
    for j in range(len(list_all)):
        image = list_all[j]
        raw_data,header = get_data_header(image)
        #print(image," ",raw_data.shape)
        data  = raw_data[0]
        x,y = data.shape
        for i in range(x):
            data[i] = data[i] - np.median(data[i][4114:4141])
        trim_data = data[124:1669,17:4113]
        out_image = patt_over_trim+image
        fits.writeto(out_image,trim_data,header,clobber=True)
        print_progress(j,length)
    return True

def create_master_bias(bias_list):
    print_spacers("Creating Master Bias")
    length = len(bias_list)
    stack = []
    for i in range(len(bias_list)):
        image = patt_over_trim+bias_list[i]
        data,header = get_data_header(image)
        stack.append(data)
        print_progress(i,length)
        #print(image)
    header["OBJECT"] = "MASTERBIAS"
    master = np.dstack(stack)
    master = np.median(master,axis=2)
    out_image = "master_bias.fits"
    fits.writeto(out_image,master,header,clobber=True)
    return out_image

def bias_correction(night):
    print_spacers("Bias Correction")
    bias,bheader = get_data_header(night.master_bias)
    list_to_correct = night.flats_wg + night.flats_ng + night.sci + night.lamp
    length = len(list_to_correct)
    for i in range(len(list_to_correct)):
        image = patt_over_trim + list_to_correct[i]
        data,header = get_data_header(image)
        bdata = data - bias
        header["COMMENT"] = "Bias Corrected Image"
        out_image = patt_bias_corr + list_to_correct[i]
        fits.writeto(out_image,bdata,header,clobber=True)
        print_progress(i,length)
    return True


def create_master_flat(flat_list):
    print_spacers("Create Master Flat For Spectroscopy")
    length = len(flat_list)
    stack = []
    for i in range(len(flat_list)):
        image = patt_bias_corr + flat_list[i]
        print(image)
        data,header = get_data_header(image)
        stack.append(data)
        print_progress(i,length)
    header["OBJECT"] = "MASTERFLAT"
    header["COMMENT"]= "Masterflat for spectroscopy"
    master = np.dstack(stack)
    master = np.median(master,axis=2)
    nmaster= normalize_flat(master)
    out_image = "master_flat_spec.fits"
    fits.writeto(out_image,nmaster,header,clobber=True)
    return out_image

def normalize_flat(master):
    x,y = master.shape
    #print(x,y,master.shape)
    black_body  = []
    x_axis      = []
    for e in range(y):
        black_body.append(np.median(master[:,e]))
        x_axis.append(e+1)
    print(len(black_body),len(x_axis))
    fitted = fit_func(x_axis,black_body,"polynomial")
    #normalizing
    for l in range(x):
        master[l,:] = master[l,:]/poly3(x_axis,*fitted)
    #fits.writeto("normaflat.fits",master,clobber=True)
    print(fitted)
    print(len(fitted))
    plt.title("Flat Normalization")
    plt.xlabel("Dispersion Direction")
    plt.ylabel("Intensity")
    plt.plot(x_axis,black_body)
    plt.plot(x_axis,poly3(x_axis,*fitted),label="Grade 3 Polynomial")
    plt.legend()
    plt.savefig("./img/flat_normalization_pol3.png",dpi=300,bbox_inches='tight')
    plt.show()
    plt.clf()
    return master

#polynonmial of second order
def poly2(x,a,b,c,x0):
    return a+b*(x-x0)+c*(x-x0)**2

def poly3(x,a,b,c,d,x0):
    return a+b*(x-x0)+c*(x-x0)**2+d*(x-x0)**3

def fit_func(x,y,func="polynomial"):
    if func == "polynomial":
        a,b,c,d,x0 = 1,1,1,1,1
        popt,pcov = curve_fit(poly3,x,y,p0=[a,b,c,d,x0])
        return popt
    elif func == "gauss":
        print("gauss")
    return True


def flat_correction(night):
    print_spacers("Master Flat Correction")
    list_to_correct = night.sci + night.lamp
    length = len(list_to_correct)
    master_flat = fits.getdata(night.master_flat)
    for i in range(len(list_to_correct)):
        image = patt_bias_corr + list_to_correct[i]
        data,header = get_data_header(image)
        cdata = data/master_flat
        header["COMMENT"] = "Flat corrected image"
        out_image = patt_flat_corr + list_to_correct[i]
        print_progress(i,length)
        fits.writeto(out_image,cdata,header,clobber=True)
    return True

#miscelaneous functions
def print_spacers(message):
    rows, columns = os.popen('stty size', 'r').read().split()
    print(rows,columns)
    if len(message)%2 == 1:
        message	= message+" "
    bar_length	= int(columns)
    bar 	= "="*bar_length
    blanks	= bar_length - 2
    blank_bar	= "="+" "*blanks +"="
    space_length	= int((blanks-len(message))/2)
    message_bar	= "=" + " " * space_length + message + " " * space_length + "="
    print(bar)
    #print(blank_bar)
    print(message_bar)
    #print(blank_bar)
    print(bar)
    
def print_progress(current,total):
    if current == total:
        sys.stdout.write("Progress {:.2%}\n".format(1.0*current/total))
    else:
        sys.stdout.write("Progress {:.2%}\r".format(1.0*current/total))
    sys.stdout.flush()
    return

if __name__ == '__main__':
    this_night = get_file_list()
    #apply overscan correction
    overscan_trim_corr(this_night.all)
    master_bias = create_master_bias(this_night.bias)
    this_night.set_master_bias(master_bias)
    #apply bias correction
    bias_correction(this_night)
    #create master flats for spectroscopy
    master_flat = create_master_flat(this_night.flats_wg)
    this_night.set_master_flat(master_flat)
    flat_correction(this_night)
    #bias_list = this_night.get_bias()
    #print(bias_list)
    #flat_1    = this_night.get_flats_wg()
    #print(flat_1)
    #flat_2      = this_night.get_flats_ng()
    #print(flat_2)
