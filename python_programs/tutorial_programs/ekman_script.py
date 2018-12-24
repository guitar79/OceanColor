#!/usr/bin/env python

from my_general_read_utilities import *
from pylab import *
import numpy as np
from math import *

#input for file name

fname = '/Users/mfc95/data/quickscat_files/QS_XWGRD3_2008036.20080401645'

#dimension of datasets and latitude/longitude ranges

ydim = 720
xdim = 1440

east = 180
west = -180
north = 90
south = -90

#Variables and constants needed for calculations

radians_per_pixel = (pi/180)*(north-south)/(ydim-1)
rads = pi/2.0 - np.arange(ydim)*radians_per_pixel
omega = 7.292 * pow(10,-5)
f = 2.0 * omega *np.sin(rads)
rho = 1028.0
rho_air = 1.2
c_d = 0.0012


#color map info

mycmap = get_cmap('spectral')
mycmap.set_bad('k')

#pixel/distance conversions

meters_per_pixel_lat = 111.2 * (1000)*(north-south)/(ydim-1)
meters_per_pixel_lon = meters_per_pixel_lat * np.cos(abs(rads))

#holding arrays established for volume transport and ekman pumping, needed in for loop

u_volume_transport = np.zeros((ydim,xdim))
v_volume_transport = np.zeros((ydim,xdim))
ekman_pumping = np.zeros((ydim,xdim))



# Read in quickscat files and break into separate components

wind_uv = read_quick_scat(fname)

u_component = wind_uv[0,:,:]
v_component = wind_uv[1,:,:]

u2_component = np.square(u_component)
v2_component = np.square(v_component)

#compute windspeed from U and V vectors

u2v2_component = u2_component + v2_component

uv_component= np.sqrt(u2v2_component)

#calculate u and v wind stress components

ustress = rho_air * (u_component/np.abs(u_component))*u2_component * c_d
vstress = rho_air * (v_component/np.abs(v_component))*v2_component * c_d

#calculate wind stress gradients for ekman pumping

[gradient_u_y, gradient_u_x] = np.gradient(ustress)
[gradient_v_y, gradient_v_x] = np.gradient(vstress)

#calculation for horizontal ekman volume transport


for i in range(720):
	if sin(rads[i]) < 0.07 and sin(rads[i]) > -0.07:
		u_volume_transport[i,:] = np.nan
		v_volume_transport[i,:] = np.nan
		ekman_pumping[i,:] = np.nan
		
		
	else:
		u_volume_transport[i,:] = vstress[i,:]/f[i]/rho
		v_volume_transport[i,:] = -(1)*(ustress[i,:]/f[i]/rho)
		ekman_pumping[i,:] = (1/(rho *f[i]))*(gradient_v_x[i,:]-gradient_u_y[i,:]) + ustress[i,:]
			
			#u_volume_transport[i,:] =  f[i]/rho
			#v_volume_transport[i,:] =  f[i]
			

figure('u stress')
imshow(ustress, cmap = mycmap, vmin = -0.5, vmax = 0.5)

figure('v stress')
imshow(vstress, cmap = mycmap, vmin = -0.5, vmax = 0.5)
		
figure('u volume transport')
imshow(u_volume_transport, cmap = mycmap,vmin = -5,vmax = 5)

figure('v volume transport')
imshow(v_volume_transport, cmap = mycmap,vmin = -5, vmax = 5)

figure('Ekman Pumping')
imshow(ekman_pumping, cmap = mycmap, vmin = -0.75, vmax = 0.75)
show()

