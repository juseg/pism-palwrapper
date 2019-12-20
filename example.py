#!/usr/bin/env python
# Copyright (c) 2018, Julien Seguinot <seguinot@vaw.baug.ethz.ch>
# GNU General Public License v3.0+ (https://www.gnu.org/licenses/gpl-3.0.txt)

"""Generate Bash scripts for example PISM run on the Alps."""

import pism_palwrapper

# make config and job scripts
pism_palwrapper.make_all(

    # configuration files
    config='alpcyc4',
    # esia2 essa1 sliding dry bedthermal lc

    # executables
    nodes=2,
    mpi_exec='',  #srun --ntasks-per-node 36',
    pism_exec='pismr',
    pism_root='.',

    # input and output files
    i_file='alps.srtm.5km.sd1.nc',
    atm_file='alps.wcnn.5km.sd1.nc',
    sd_file='alps.erai.5km.sd1.nc',
    dt_file='epica.3222.1220.sd1.nc',
    out_dir='output/example',

    # vertical grid
    lbz=3000.0,
    lz=5000.0,
    mbz=13,
    mz=51,

    # time steps
    ys=-120000,
    ye=0,
    ychain=20000,
    yextra=100,
    yts=10,
    submit=False,

)
