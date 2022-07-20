# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 11:05:09 2022

@author: msnow
"""

import numpy as np
import astropy
from astropy.io import fits
import glob
import sunpy.map
import matplotlib.pyplot as plt
#from sunpy.net import Fido, attrs as a

#result=Fido.search(a.Instrument.IRIS)


iris_dir='C:/Users/Thulisile Dlamaka/Downloads'


flist=glob.glob(iris_dir+'*_HMI.fits')
with fits.open(flist[0]) as hdul_hmi:
    #hdul_hmi.info()
    hmi_data=hdul_hmi[0].data

hmi_center_x=hdul_hmi[0].header['crpix1']
hmi_center_y=hdul_hmi[0].header['crpix2']
hmi_stripe=hmi_data[:,int(hmi_center_y)]  #value in header is a float
hmi_dx=hdul_hmi[0].header['cdelt1']
hmi_x_vector=(np.arange(4096)-hmi_center_x)*hmi_dx

# with fits.open(iris_dir+'outfile.fits') as hdul:
#     hdul.info()
#     iris_data=hdul[0].data
#     hmi_data=hdul[3].data
    
# print(repr(hdul[3].header))
    
#aia=sunpy.map.Map(flist[0])
#d,h = astropy.io.fits.read(flist[0])[0]
#c=sunpy.map.all_coordinates_from_map(aia)
#print(c[0:5,0:5])

#aia.peek(data)
# fig=plt.figure(1,figsize=[10,7])
# ax=plt.subplot(111,projection=aia)
# aia.plot()
# aia.draw_limb()
# aia.draw_grid()
# plt.colorbar()
# plt.show


flist_iris=glob.glob(iris_dir+'IRIS.fits.gz')
with fits.open(flist_iris[0]) as hdul2:
    #hdul2.info()
    iris_data=hdul2[0].data
    iris_mean_profile=hdul2[1].data
    

    
w0=hdul2[0].header['crval3'] #center wavelength
dw=hdul2[0].header['cdelt3'] #delta lambda
w=(np.arange(101)-50)*dw+w0
w=w/10. #convert to nm

center_pixel_x=hdul2[0].header['crpix2']
center_pixel_y=hdul2[0].header['crpix1']

iris_total=np.zeros((hdul2[0].header['naxis1'],hdul2[0].header['naxis2']))
irisx=hdul2[0].header['naxis1']
irisy=hdul2[0].header['naxis2']
for x in np.arange(irisx):
    for y in np.arange(irisy):
        iris_total[x,y]=sum(iris_data[:,y,x])
#for idx in np.ndenumerate(iris_total):
 #   iris_total[idx]=sum(iris_data[:,idx])


#make some plots!
plt.figure(2)
plt.plot(w,iris_data[:,center_pixel_x,center_pixel_y])
plt.xlabel('Wavelength [nm]')
plt.ylabel('Counts')
plt.plot(w,iris_mean_profile[:],color='red')

plt.figure(3)
plt.plot(hmi_x_vector,abs(hmi_stripe))
plt.xlabel('Distance from disk center in arc seconds')
plt.ylabel('Magnetic field strength')
plt.title('HMI disk center profile')

plt.figure(4)
plt.contour(hmi_data,levels=[0,10,50,100])
#plt.title=('HMI Contour Map')

plt.figure(5)
plt.plot(iris_total[:,3000])

plt.figure(6)
plt.contour(iris_total,levels=[5000,10000,15000,20000])

plt.show()






