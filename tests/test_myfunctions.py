import map_res
import numpy as np

e4d_inp_f = 'emesh7t'
xnods = np.arange(10)
ynods = np.arange(10)
znods = np.arange(10)

map_res.mesh_interp(xnods,ynods,znods,e4d_inp_f) # create'.bin' executable for interpolation


nv = (len(xnods)-1)*(len(ynods)-1)*(len(znods)-1)
fcr = np.zeros((nv,1))+0.4 
porosity = np.zeros((nv,1))+0.3 
temperature = np.zeros((nv,1))+25.0
satr = np.zeros((nv,1))+ 1.0 # fully saturated

sigfile='baseline.sig'
petfile='wax_smit.txt'
time = 0.0
mapfile = e4d_inp_f + '_map.bin'
map_res.map_waxsmit(fcr,satr,porosity,temperature,sigfile,petfile,mapfile,time)       
