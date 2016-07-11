#!/usr/bin/python3.4
"""
Simon Torres 2016-06-28
pipeline for GOODMAN spectra reduction.
"""



import sys
import os
import glob
import numpy as np
import re
import time
import astropy.stats as asst
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')


class night:
    """
    The night class stores the data relative to single observing night
    therefore this software works on a per-night basis
    """

    def __init__(self,date):
        self.all     = []
        self.date    = date
        self.sci     = []
        self.lamp    = []

    def add_sci(self,insci):
        """Adds science object to list"""
        self.sci.append(insci)
        self.all.append(insci)
    
    def add_lamp(self,inlamp):
        """Adds lamp objects to list"""
        self.lamp.append(inlamp)
        self.all.append(inlamp)

    
    
class sci_obj:
    """class that defines a science object atributes

    Sci objects, for science, are the main targets which have lamps for
     calibration. Their atritutes are: the name, the file name, the
     observation time, right ascension and declination. The same information
     can be obtained for the lamps but they are treated as lists
     The list are created in a way that allows the use of elements' index
     as the correlator between atributes of the lamp.

     Attributes:
         name (str): science object name
         file_name (str): file name
         obs_time (str): observing time in the format yyyy-mm-ddThh:mm:ss.ss for instance 2016-03-20T23:54:15.96
         ra (float): right ascension in degrees
         dec (float): declination in degrees
         lamp_count (int): lamps count
         lamp_file (list): every element is a string with the file name of the lamp
         lamp_type (list): every element is a string with the OBJECT value of the lamp i.e Cu, HgAr, etc
         lamp_ra (list): every element is a float with lamp's right ascension in degrees
         lamp_dec (list): every element is a float with lamp's declination in degrees

    """
    
    def __init__(self,name,file_name,obs_time,ra,dec):
        self.name       = name
        self.file_name  = file_name
        self.obs_time   = obs_time
        self.ra         = ra
        self.dec        = dec
        self.lamp_count = 0
        self.lamp_file  = []
        self.lamp_type  = []
        self.lamp_ra    = []
        self.lamp_dec   = []
    
    def add_lamp(self,new_lamp,new_type,ra,dec):
        """Adds a lamp to the science object

        Args:
            new_lamp (str): new lamp file name
            new_type (str): new lamp type
            ra (float): right ascension in degrees
            dec (float): declination in degrees

        """
        self.lamp_file.append(new_lamp)
        self.lamp_type.append(new_type)
        self.lamp_ra.append(ra)   
        self.lamp_dec.append(dec)
        self.lamp_count = int(len(self.lamp_file))

    def print_all(self):
        """Prints all the relevant atributes of the object

        Note:
            this method is mainly used for development purposes
        """
        print("Name: ",self.name)
        print("File: ",self.file_name) 
        print("Obs-T: ",self.obs_time)
        if self.lamp_count > 0:
            print("Lamp N: ",self.lamp_count)
            for i in range(self.lamp_count):
                print("Lamp %s: "%(i+1),self.lamp_file[i])
                print("Type %s: "%(i+1),self.lamp_type[i])
        
        
    
    

def get_file_list():
    """Gets flatfielded science and lamp targets

    Args:
        No arguments required.

    Returns:
        Night class object which contains separated lists
        for science and lamp targets.

    In case an IOError exception is raised prints an
    error message and exit the program.
    """

    lista = sorted(glob.glob("fc_0*.fits")) #will need work
    try:
        header0 = fits.getheader(lista[0])
        this_night = night(header0["DATE"])
        for i in lista:
            #print("./" +i)
            header  = fits.getheader(i)
            obstype = header["OBSTYPE"]
            if obstype == "OBJECT":
                this_night.add_sci(i)
            elif obstype == "COMP":
                this_night.add_lamp(i)
        return this_night
    except IOError as err:
        if str(err) == "Empty or corrupt FITS file":
            sys.exit("Raised an IOError as ",err)
        else:
            print(err)
            sys.exit("Please correct the errors and try again")
        
#GET DATA
def get_data_header(file_name):
    """Get data and header of a fits image

    Args:
        file_name (str): Full or relative path to the file to open

    Returns:
        data and header
    """
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
    
    
def convert_time(in_time):
    return time.mktime(time.strptime(in_time,"%Y-%m-%dT%H:%M:%S.%f"))

def ra_dec_to_deg(ra,dec):
    ra = ra.split(":")
    dec = dec.split(":")
    #RIGHT ASCENTION conversion
    ra_deg = (float(ra[0])+(float(ra[1])+(float(ra[2])/60.))/60.)*(360./24.)
    #DECLINATION conversion
    sign = float(dec[0])/abs(float(dec[0]))
    dec_deg = sign * (abs(float(dec[0]))+(float(dec[1])+(float(dec[2])/60.))/60.)
    return ra_deg,dec_deg

    
def organize_extraction(night):
    sci_list    =  night.sci
    lamp_list   =  night.lamp
    #loops through science objects
    for i in range(len(sci_list)):
        """
        Captures science object metadata
        RA and DEC are stored in degrees for comparison
        TIME is stored in seconds since epoch
        """
        header      = fits.getheader(sci_list[i])
        name        = header["OBJECT"]
        exptime     = header["EXPTIME"]
        file_name   = sci_list[i]
        raw_time    = header["DATE-OBS"]
        obs_time    = convert_time(raw_time)
        ra,dec      = ra_dec_to_deg(header["RA"],header["DEC"])
        
        #print(sci_list[i],obs_time)
        #create sci_obj object
        sci = sci_obj(name,file_name,raw_time,ra,dec)
        #loops through lamps
        time_diff   = []
        lamp_file   = []
        lamp_type   = []
        lamp_ra     = []
        lamp_dec    = []
        #lamp_dist   = []
        for e in range(len(lamp_list)):
            lamp        = lamp_list[e]
            lamp_header = fits.getheader(lamp)
            lamp_time   = lamp_header["DATE-OBS"]
            lamp_obj    = lamp_header["OBJECT"]
            lamp_ra_val,lamp_dec_val = ra_dec_to_deg(lamp_header["RA"],lamp_header["DEC"])
            lamp_dist   = np.sqrt((ra-lamp_ra_val)**2 + (dec-lamp_dec_val)**2)
            time_diff.append(abs(obs_time-convert_time(lamp_time)))
            lamp_file.append(lamp)
            lamp_type.append(lamp_obj)
            lamp_ra.append(lamp_ra_val)
            lamp_dec.append(lamp_dec_val)
            if lamp_dist < 1e-3:
                if time_diff[e] < exptime+900:
                    #print(time_diff[e],exptime+300)
                    #print("Add lamp: ",lamp)
                    sci.add_lamp(lamp,lamp_obj,lamp_ra_val,lamp_dec_val)
                else:
                    #print(time_diff[e],exptime,exptime+300,lamp,lamp_dist)
                    print("Lamp-to-object difference is out of range ")
        if len(sci.lamp_file) == len(sci.lamp_type):
            #print("\n\n")
            if len(sci.lamp_file) < 1:
                print("Warning! Object %s has no calibration lamps"%sci.name)
            else:
                sci.print_all()
                get_spectrum(sci)
                print("\n\n")
        else:
            sys.exit("Mismatch of lamps")
    return True

def  get_spectrum(sci_object):
    print_spacers(sci_object.name)
    data,header = get_data_header(sci_object.file_name)
    x,y = data.shape
    sample_point= np.linspace(0,y,100,dtype=int)
    sample_point= sample_point[1:-1]
    y_max   = []
    x_max   = []
    for i in sample_point:
        """   Need to identify secondary objects   """
        y_max.append(np.argmax(data[:,i]))
        x_max.append(i)
        #print(np.argmax(data[:,i]))
    #sigmaclip to remove outliers
    y_max = asst.sigma_clip(y_max,sig=3,cenfunc=np.median)
    """fit gaussian to slit profile """
    gauss_check = np.linspace(0,len(x_max),10,dtype=int)
    gauss_x = []
    gauss_y = []
    for e in range(len(x_max)):
        if e in gauss_check:
            gauss_x.append(x_max[e])
            gauss_y.append(y_max[e])
            #print(x_max[e],y_max[e])
    gauss_sigma =[]
    for g in range(len(gauss_x)):
        gy      = data[:,gauss_x[g]]
        gx      = range(len(gy))
        mean    = gauss_y[g]
        offset  = mean
        sigma   = 4
        popt,pcov = curve_fit(gauss,gx,gy,p0=[1.,mean,sigma,offset])
        gauss_sigma.append(popt[2])
    gauss_sigma_length = len(gauss_sigma)
    gauss_sigma = asst.sigma_clip(gauss_sigma,sig=2,cenfunc=np.median)
    sigma       = np.mean(gauss_sigma)
    for g in range(len(gauss_x)):
        ngy     = data[int(y_max[g]-4*sigma)-1:int(y_max[g]+4*sigma),gauss_x[g]]
        ngx     = range(len(ngy))
        #interpolate
        igx     = np.linspace(ngx[0],ngx[-1],len(ngx)*100)
        tck     = interpolate.splrep(ngx,ngy,s=0)
        igy     = interpolate.splev(igx,tck,der=0)
        
        mean    = gauss_y[g]
        offset  = mean
        popt,pcov   = curve_fit(gauss,igx,igy,p0=[1.,mean,sigma,offset])
        print(popt)
        plt.title("Column: %s"%(gauss_x[g]))
        plt.xlabel("spatial Direction")
        plt.ylabel("Intensity")
        plt.plot(gx,gy)
        plt.plot(gx,gauss(gx,*popt),label="Gaussian Fit")
        plt.legend()
        plt.savefig("./img/gauss_fit_"+sci_object.name+"L"+str(gauss_x[g])+"-2.png",dpi=300,bbox_inches='tight')
        #plt.show()
        plt.clf()
        
    #figure parameters
    plt.title(sci_object.name)
    plt.xlabel("Dispersion Direction")
    plt.ylabel("Spatial Direction")
    plt.imshow(data,cmap='gray',clim=(5,150))
    plt.plot(x_max,y_max,lw=0.5,color='r')
    plt.savefig("./img/"+sci_object.name+'-trace.png',dpi=300,bbox_inches='tight')
    #plt.show()
    plt.clf()
    #end of figure
    print(np.std(y_max))
    print(y_max)
    print("median ",np.median(data))
    return True

#def normalize_flat(master):
    #x,y = master.shape
    ##print(x,y,master.shape)
    #black_body  = []
    #x_axis      = []
    #for e in range(y):
        #black_body.append(np.median(master[:,e]))
        #x_axis.append(e+1)
    #print(len(black_body),len(x_axis))
    #fitted = fit_func(x_axis,black_body,"polynomial")
    ##normalizing
    #for l in range(x):
        #master[l,:] = master[l,:]/poly3(x_axis,*fitted)
    ##fits.writeto("normaflat.fits",master,clobber=True)
    #print(fitted)
    #print(len(fitted))
    #plt.plot(x_axis,black_body)
    #plt.plot(x_axis,poly3(x_axis,*fitted),label="fit")
    #plt.legend()
    #plt.show()
    #return master

#polynonmial of second order
def poly2(x,a,b,c,x0):
    return a+b*(x-x0)+c*(x-x0)**2

def poly3(x,a,b,c,d,x0):
    return a+b*(x-x0)+c*(x-x0)**2+d*(x-x0)**3

def gauss(x,a,x0,sigma,c):
    return c+a*np.exp(-(x-x0)**2/(2*sigma**2))

def fit_func(x,y,func="polynomial"):
    if func == "polynomial":
        a,b,c,d,x0 = 1,1,1,1,1
        popt,pcov = curve_fit(poly3,x,y,p0=[a,b,c,d,x0])
        return popt
    elif func == "gauss":
        print("gauss")
    return True




#miscelaneous functions
def print_spacers(message):
    rows, columns = os.popen('stty size', 'r').read().split()
    if len(message)%2 == 1 and int(columns)%2 != 1:
        message	= message+" "
    bar_length	= int(columns)
    bar 	= "="*bar_length
    blanks	= bar_length - 2
    blank_bar	= "="+" "*blanks +"="
    space_length	= int((blanks-len(message))/2)
    message_bar	= "=" + " " * space_length + message + " " * space_length + "="
    print(bar)

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
    #organize observations
    organize_extraction(this_night)
    print(this_night.sci)
    print(this_night.lamp)
    #bias_list = this_night.get_bias()
    #print(bias_list)
    #flat_1    = this_night.get_flats_wg()
    #print(flat_1)
    #flat_2      = this_night.get_flats_ng()
    #print(flat_2)
