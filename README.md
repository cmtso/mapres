# map_res
A lightweight efficient python library to interpolate and map output from a 3D flow and transport model (rectangular mesh) to a ERT mesh and perform petrophysical transform to obtain resistivities. It is useful for coupled hydrogeophysical simulations to link flow and transport code to ERT code. An example application is mapping from SEAWAT to RESIPY, which I include in a accompanying Jupyter notebook.

The library includes two FORTRAN subroutines written by Tim Johnson for [PFLOTRAN-E4D](https://doi.org/10.1016/j.cageo.2016.09.006) (Johnson et al. 2017). It was part of the PFLOTRAN distribution but has now discontinued. He has kindly provided the FORTRAN subroutines and some other python scripts and has given me permission to package it as a python library here. Please cite the original PFLOTRAN-E4D paper and find more details of the method there.

FORTRAN codes can be called in python using  `f2py` in `numpy`. In its original form, the user needs to compile the FORTRAN subroutines on their own before they can use it. I can packaged packaged it as a python library for convenience (so that the f2py installation steps are automated and there is no need to track the path of the scripts).

**You must cite this work using the instruction below.**



Installation 
-------------------
- Your machine needs to contain a fortran compiler, e.g. gfortran. On Linux debian/ubuntu, you can obtain it by running:

```sh
sudo apt update
sudo apt-get install gfortran
```
- Download or clone this repo. Go to the folder.
- You can install in the Terminal or Anaconda prompt:
```sh
python -m pip install .
```


Test the installation
---------------------
After installation, run
```sh
cd tests
python test_myfunctions.py
```
If terminated normally, two files *emesh7t.map* and *0.00000000.sig* are generated.

Examining the [test file](tests/other_file.md) will help you understand how to use this library.

What it does and how it works?
------------------------------------------
It contains two functions:
- `map_res.somef90.mesh_interp(xnods,ynods,znods,e4d_inp_f)` create a '.bin' binary executable for interpolation. It determines to weights for interpolation and is the time-consuming part. It only needs to run once even if you need to map multiple timesteps, as long as the two meshes stayed unchanged.
- `map_res.somef90.map_waxsmit(...)` maps variables (e.g. fluid conductivity and saturation) very efficiently between the two meshes once the binary executable is created and return **electrical conductivity (EC, unit=S/m)** for each cell in the ERT mesh via petrophysical transform (Archie/Waxman Smit).


API reference and input files formats
------------------------------------------

- `xnods,ynods,znods`: x,y,z coordinates of the flow and transport model in the x,y,z direction respectively. Note the only simple rectangular grid is supported and they should be strictly increasing.
- `e4d_inp_f`: The ERT mesh file prefixes (in tetgen or E4D format). The program will read a series of mesh files `.1.node`(node list), `.1.ele` (cell list, what node constitutes each element), `.1.face`,`.1.neigh`,`.1.trn` (translation, [0 0 0] if unchanged)

Note that [RESIPY](https://gitlab.com/hkex/resipy) has support to generate and read E4D file formats. Consult the E4D documentation if needed.

- `fcr`: an 1D array of fluid conductivity in each cell of the flow and transport model
- `satr`: an 1D array of saturation in each cell of the flow and transport model
- `temperature`: an 1D array of fluid conductivity in each cell of the flow and transport model
- `sigfile`: file for baseline electrical conductivity for each cell of the ERT model. ERT cells not interpolated will use this value.
- `petfile`: file for petrophysical parameters for each cell of the ERT model (see below)
- `mapfile`: name of the binary mapping file. Just use `e4d_inp_f+'_map.bin'`
- `time`: the timestep to tag the output file name. 

The output file `<time>.00000000.sig` contains the EC values for each cell of the ERT mesh, which can then be loaded to a ERT code to run as a forward model.

Ususally, flow and transport variables are outputted as 3D arrays of size . Assuming `fcond` and `sat` are variables in 3D array format, they can be converted to 1D arrays using
```python
import numpy as np
nv = (xnods.shape[0]-1)*(ynods.shape[0]-1)*(znods.shape[0]-1) # number of cells for flow and transport grid
fcr = np.reshape(fcond,(nv,1),order='F')
satr = np.reshape(sat,(nv,1),order='F')
```


Petrophysical transform
-----------------------
The program reads Waxman-Smit parameters for each cell. The columns are 1) a, 2) B, 3) Qv, 4) c (cementation factor), 5) m (saturation exponent), 6) t(temperature correction factor).

To use Archie transform, just set $a=1, B=0, Qv=0$.

For other etrophysical transforms, you can modify the code accordingly.


Other applications of this mapping script
------------------------------------------

- eSTOMP to E4D (Rockhold et al. 2020) https://doi.org/10.1007/s10040-020-02167-1
- PFLOTRAN to E4D (Tso et al. 2019) https://doi.org/10.1007/s10040-020-02167-1


You must cite:
------------------
If you use this library for your work, please cite [this paper](https://doi.org/10.1016/j.cageo.2016.09.006) as:

    Timothy C. Johnson, Glenn E. Hammond, Xingyuan Chen,
    PFLOTRAN-E4D: A parallel open source PFLOTRAN module for simulating time-lapse electrical resistivity data,
    Computers & Geosciences,
    Volume 99,
    2017,
    Pages 72-80,
    https://doi.org/10.1016/j.cageo.2016.09.006.

BibTex code:
```latex
@article{JOHNSON201772,
    title = {{PFLOTRAN-E4D: A parallel open source PFLOTRAN module for simulating time-lapse electrical resistivity data}},
    journal = {{Computers & Geosciences}},
    volume = {99},
    pages = {72 - 80},
    year = {2017},
    doi = {https://doi.org/10.1016/j.cageo.2016.09.006},
    author = {imothy C. Johnson and Glenn E. Hammond and Xingyuan Chen},
    keywords = {Hydrogeophysics, Time-lapse geophysics, Electrical resistivity tomography, Groundwater, Simulation, Multi-physics, Parallel, Open-source}
}
```

