#!/usr/bin/env python2

"""A PISM wrapper for paleo-jobs."""

import os
import subprocess

# global system settings
mpi_exec = 'aprun -B'
pism_version = '0.7.2'
pism_arch = 'dora-intel'
pism_exec = os.path.expanduser('~/software/opt/pism-%s/%s/bin/pismr'
                               % (pism_version, pism_arch))
pism_root = os.path.expanduser('~/pism')


# job script template
# FIXME: add slurm magic comments
# FIXME: make topg_to_phi an option
template = '''#!/bin/bash
{mpi_exec} {pism_exec} \\
    -i {boot_path} -bootstrap \\
    -config_override config.nc -topg_to_phi 15,45,0.0,200.0 \\
    -Mx {mx} -My {my} -Mz {mz} -Mbz {mbz} -Lz 5000 -Lbz 3000 \\
    -ys {ys} -ye {ye} -z_spacing equal \\
    -atmosphere given,lapse_rate,delta_T -temp_lapse_rate 6.0 \\
    -atmosphere_given_file {atm_path} \\
    -atmosphere_given_period 1 -timestep_hit_multiples 1 \\
    -atmosphere_lapse_rate_file {atm_path} \\
    -atmosphere_delta_T_file {dt_path} \\
    -surface pdd -pdd_sd_period 1 \\
    -pdd_sd_file {sd_path} \\
    -ocean pik,delta_SL \\
    -ocean_delta_SL_file {dsl_path} \\
    -ts_file {prefix}-ts.nc -ts_times {yts} \\
    -extra_file {prefix}-extra.nc -extra_times {yextra} \\
    -extra_vars bmelt,climatic_mass_balance,cbase,csurf,lat,lon,mask,rank,\\
tauc,taud_mag,tempicethk_basal,temppabase,tempsurf,thk,topg,usurf,\\
velbase,velbase_mag,velsurf,velsurf_mag \\
    -o {prefix}.nc > {prefix}.log 2> {prefix}.err
'''


# set region extents in km
# FIXME: extract this info from boot files instead
extents = {
    # North America
    'cordillera'    : (1500, 3000),
    'olympic'       : ( 192,  192),
    'puget'         : ( 600,  750),
    # Europe
    'alps'          : ( 900,  600),
    # Asia
    'haizishan'     : (  75,  125),
    'kodar'         : ( 360,  360),
    'stanovoy'      : ( 600,  400),
    'transbaikalia' : (1000,  800)}


def make_config(config, out_dir=None):
    """Create configuration file and return its path."""
    # FIXME: read config files as dicts and export to netCDF
    # FIXME: allow to concatenate multiple files 
    cdl_path = '%s/config/%s.cdl' % (pism_root, config)
    nc_path = os.path.join(out_dir, 'config.nc')
    cmd = ['ncgen', cdl_path, '-o', nc_path]
    subprocess.call(cmd)
    return nc_path


def make_jobscript(reg, res, boot_file, atm_file, sd_file, dt_file, dsl_file,
                   mz=51, mbz=31, ys=0.0, ye=1000.0, yts=10, yextra=100,
                   out_dir=None):
    """Create job script and return its path."""

    # set region extents in km
    x, y = extents[reg]

    # compute number of grid points
    mx = x*1000/res + 1
    my = y*1000/res + 1
    if atm_file.endswith('.cr.nc'):
        mx -= 1
        my -= 1

    # parse paths
    boot_path = os.path.join(pism_root, 'input', 'boot', boot_file)
    atm_path = os.path.join(pism_root, 'input', 'atm', atm_file)
    sd_path = os.path.join(pism_root, 'input', 'sd', sd_file)
    dt_path = os.path.join(pism_root, 'input', 'dt', dt_file)
    dsl_path = os.path.join(pism_root, 'input', 'dsl', dsl_file)

    # prefix for output files
    prefix = 'y%07d' % (ye-ys)

    # format script
    script = template.format(mpi_exec=mpi_exec, pism_exec=pism_exec,
                             **locals())

    # write script to file
    script_path = os.path.join(out_dir, prefix + '.sh')
    with open(script_path, 'w') as f:
        f.write(script)

    # return path to job script
    return script_path


def make_all(reg, res, boot_file, atm_file, sd_file, dt_file, dsl_file, config,
             out_dir, **kwargs):
    """Create new directory, job script and config file."""

    # make new directory
    # FIXME: recursive mkdir
    # FIXME: better handle existing dir case
    os.mkdir(out_dir)

    # make config file
    c_path = make_config(config, out_dir=out_dir)

    # make job script
    j_path = make_jobscript(reg, res,
                            boot_file, atm_file, sd_file, dt_file, dsl_file,
                            out_dir=out_dir, **kwargs)

    # print path to new jobscript
    print j_path
