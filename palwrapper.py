#!/usr/bin/env python2

"""A PISM wrapper for paleo-jobs."""

import os
import subprocess
from netCDF4 import Dataset

# global system settings
mpi_exec = 'aprun -B'
pism_version = '0.7.2'
pism_arch = 'dora-intel'
pism_exec = os.path.expanduser('~/software/opt/pism-%s/%s/bin/pismr'
                               % (pism_version, pism_arch))
pism_root = os.path.expanduser('~/pism')


# job script template
template = '''#!/bin/bash
#
#SBATCH --job-name={prefix}
#SBATCH --nodes={nodes}
#SBATCH --time={time}
#SBATCH --output={prefix}.log
#SBATCH --error={prefix}.err

{mpi_exec} {pism_exec} \\
    -i {boot_path} {boot_args} \\
    -o {prefix}.nc \\
    -ys {ys} -ye {ye}\\
    -config_override config.nc \\
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
velbase,velbase_mag,velsurf,velsurf_mag

'''

# FIXME: make topg_to_phi an option
boot_args_template = '''\\
    -bootstrap -Mx {mx} -My {my} -Mz {mz} -Mbz {mbz} -Lz 5000 -Lbz 3000 \\
    -z_spacing equal -topg_to_phi 15,45,0.0,200.0 '''


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
    # FIXME: allow to concatenate multiple files 

    # input and output file paths
    txt_path = '%s/config/%s.txt' % (pism_root, config)
    nc_path = os.path.join(out_dir, 'config.nc')

    # initialize netCDF dataset
    nc = Dataset(nc_path, 'w')
    var = nc.createVariable('pism_overrides', 'i1')

    # fill in pism overrides
    with open(txt_path) as f:
        for line in f:
            k, v = line.rstrip().split(None, 1)
            v = v.strip('"')
            if not k.startswith('//'):
                var.setncattr(k, v)

    # close and return path to output file
    nc.close()
    return nc_path


def make_jobscript(reg, res, boot_file, atm_file, sd_file, dt_file, dsl_file,
                   mz=51, mbz=31, ys=0.0, ye=1000.0, yts=10, yextra=100,
                   nodes=1, time='24:00:00', out_dir=None, prefix='run',
                   bootstrap=True):
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
    atm_path = os.path.join(pism_root, 'input', 'atm', atm_file)
    sd_path = os.path.join(pism_root, 'input', 'sd', sd_file)
    dt_path = os.path.join(pism_root, 'input', 'dt', dt_file)
    dsl_path = os.path.join(pism_root, 'input', 'dsl', dsl_file)

    # bootstrapping arguments
    if bootstrap == True:
        boot_path = os.path.join(pism_root, 'input', 'boot', boot_file)
        boot_args = boot_args_template.format(**locals())
    else:
        boot_path = boot_file
        boot_args = ''

    # format script
    script = template.format(mpi_exec=mpi_exec, pism_exec=pism_exec,
                             **locals())

    # write script to file
    script_path = os.path.join(out_dir, prefix + '.sh')
    with open(script_path, 'w') as f:
        f.write(script)

    # return path to job script
    return script_path


def make_chain(reg, res, boot_file, atm_file, sd_file, dt_file, dsl_file,
               **kwargs):
    """Create several job scripts to run as a chain."""

    # pop relevant keyword arguments
    ys = kwargs.pop('ys', 0.0)
    ye = kwargs.pop('ye', 1000.0)
    ychain = kwargs.pop('ychain', None)

    # set ychain to match run duration if invalid
    if ychain is None or ychain <= 0.0 or ychain > (ye-ys):
        ychain = ye - ys

    # create the first jobscript
    boot_job_name = 'y%07d' % (ychain)
    boot_job_path = make_jobscript(reg, res,
                                   boot_file, atm_file, sd_file, dt_file, dsl_file,
                                   ys=ys, ye=ys+ychain, prefix=boot_job_name,
                                   bootstrap=True, **kwargs)
    job_path_list = [boot_job_path]

    # create the next jobscripts if necessary
    i_file = boot_job_name + '.nc'
    if ychain < (ye-ys):
        for y in range(ys+ychain, ye, ychain):
            job_name = 'y%07d' % (ychain+y-ys)
            job_path = make_jobscript(reg, res,
                                      i_file, atm_file, sd_file, dt_file, dsl_file,
                                      ys=y, ye=y+ychain, prefix=job_name,
                                      bootstrap=False, **kwargs)
            job_path_list.append(job_path)
            i_file = job_name + '.nc'

    # return path to first job script
    return job_path_list


def submit_job(job_path, depends=None):
    """Submit a job script and return job ID."""

    # tell sbatch to work in dir containing script
    job_dir = os.path.dirname(job_path)
    cmd = ['sbatch', '--workdir='+job_dir]

    # add dependency if any
    if depends is not None:
        cmd.append('--dependency=afterok:'+depends)

    # run sbatch command
    cmd.append(job_path)
    out = subprocess.check_output(cmd)

    # return job id
    job_id = out.rstrip('\n').split(' ')[-1]
    return job_id


def submit_chain(job_path_list):
    """Submit a list of job scripts as a chain and return job IDs."""

    # run first job
    job_id = submit_job(job_path_list[0])
    job_id_list = [job_id]

    # run other jobs if any
    for job_path in job_path_list[1:]:
        job_id = submit_job(job_path, depends=job_id)
        job_id_list.append(job_id)

    # return list of job ids
    return job_id_list


def make_all(reg, res, boot_file, atm_file, sd_file, dt_file, dsl_file, config,
             out_dir, submit=True, **kwargs):
    """Create new directory, job script and config file."""

    # make new directory or break if existing
    try:
        os.makedirs(out_dir)
    except OSError:
        print "Directory %s exists, skipping it." % out_dir
        return 2

    # make config file
    c_path = make_config(config, out_dir=out_dir)

    # make job script chain
    j_list = make_chain(reg, res,
                        boot_file, atm_file, sd_file, dt_file, dsl_file,
                        out_dir=out_dir, **kwargs)

    # submit job chain and print job ids
    if submit == True:
        j_list = submit_chain(j_list)
        print 'Submitted jobs: ' + ' '.join(j_list)

    # or print list of scripts
    else:
        print 'Create scripts:\n' + '\n'.join(j_list)

    # no error, return 0
    return 0
