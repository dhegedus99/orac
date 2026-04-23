"""Routines to run an ORAC component."""
import os, re, datetime
import pyorac.arguments as oracarg
import pyorac.definitions as defin
from pyorac.colour_print import colour_print

from copy import deepcopy
from collections import OrderedDict
from glob import glob
from pyorac import defaults
from pyorac.util import call_exe, read_orac_library_file, build_orac_library_path
from subprocess import check_call, check_output, CalledProcessError
import tempfile

CLOBBER = OrderedDict([
    ('pre', 5),
    ('main', 4),
    ('post', 3),
    ('flux', 2),
    ('sisem', 1)
])

def process_latest_ecmwf_files(tartime, rawecm, proecm, ecmwf_in, ecmwf_out, ecsdir, 
                               args, tag='ecmwf', dependency=None):
    from pyorac import convertgrb2nc as grb2nc
    rawecm1 = [x for x in rawecm if re.match(r'\w{11,11}'+tartime+'\w',x)]
    proecm1 = [x for x in proecm if re.match(r'\w{11,11}'+tartime+'\w*',x) and x.endswith('.nc')]
    process_ecmwf = ''
    delete_ecmwf = ''
    # Make sure we're not looking at data produced at the time step 
    # for which it is forecasting (i.e. the first and second time
    # stamps are the same). Such data aren't the same as other files
    if abs(int(max(rawecm1)[11:19]) - int(max(rawecm1)[3:11])) < 100:
        rawecm1 = rawecm1[0:len(rawecm1)-1]
    if len(rawecm1) > 0 and len(proecm1) > 0 and \
            (max(rawecm1))[3:11] > (max(proecm1))[3:11]:                    
        process_ecmwf=max(rawecm1)
        delete_ecmwf=max(proecm1)
    elif len(rawecm1) > 0 and len(proecm1) == 0:
        process_ecmwf=max(rawecm1)
    if len(process_ecmwf) > 0:
        if args.batch:
            job_name = args.File.job_name(args.revision, tag)
            from pyorac.local_defaults import LOG_DIR
            import string, random
            log_path = os.path.join(args.out_dir, LOG_DIR)
            uid = ''.join([random.choice(string.ascii_letters+string.digits) 
                   for n in range(6)])
            values = {'job_name': job_name,
                    'log_file': os.path.join(log_path, job_name + uid +'.log'),
                    'err_file': os.path.join(log_path, job_name + uid +'.err'),
                    'account': defaults.BATCH_VALUES['account'],
                    'qos' : defaults.BATCH_VALUES['qos'],
                    'queue': defaults.BATCH_VALUES['queue']}
            if 'ram' in defaults.BATCH_VALUES:
                values['ram'] = defaults.BATCH_VALUES['ram'][0]
            else:
                values['ram'] = args.ram[0]
            if 'duration' in defaults.BATCH_VALUES:
                values['duration'] = defaults.BATCH_VALUES['duration'][0]
            else:
                values['duration'] = args.dur[0]
            if dependency is not None:
                values['depend'] = dependency
            exe = args.orac_dir+'/tools/pyorac/convertgrb2nc.py'
            cmd = 'python ' + exe + ' ' \
                    + ecmwf_in + ecsdir + ' ' + ecmwf_out + ecsdir + ' ' + process_ecmwf
            # Write temporary script to call executable
            (gd, script_file) = tempfile.mkstemp('.sh', os.path.basename(exe)+'.',
                                             args.out_dir, True)
            g = os.fdopen(gd, "w")
            g.write("#!/bin/bash\n")
            # Define processing environment
            libs = read_orac_library_file(args.orac_lib)
            g.write("export LD_LIBRARY_PATH=" +
                              build_orac_library_path(libs) + "\n")
            g.write("export OPENBLAS_NUM_THREADS=1\n")
            defaults.BATCH.add_openmp_to_script(g)
            
            g.write(cmd+"\n")
            g.write("rm -f "+script_file+"\n")
            g.close()
            os.chmod(script_file, 0o700)
            cmd = defaults.BATCH.list_batch(values, exe=script_file)
            if args.verbose or args.script_verbose:
                    colour_print(' '.join(cmd), defin.COLOURING['header'])
            out = check_output(cmd, universal_newlines=True)
    
            # Parse job ID # and return it to the caller
            jid = defaults.BATCH.parse_out(out, 'ID')
        else:
            status = grb2nc.proc_dir(ecmwf_in+ecsdir, ecmwf_out+ecsdir, process_ecmwf)
            if status != 0:
                print('Warning: Convert_ECM_GRB2NC.py encountered an ', \
                    'error with one or more files on date: ',ecsdir)
            else:
                jid=True
        if len(delete_ecmwf) > 0:
            print('to delete', delete_ecmwf)
        ecm_out = ecmwf_out+ecsdir + '/' + process_ecmwf +'.nc'
    else:
        jid = None
        ecm_out = None
    return jid, ecm_out

def pre_process_ecmwf_grib(yr, mth, day, hr, ecmwf_out, ecmwf_in_1, args):
    # If we're on, or after, the final ECMWF timestep of the day, we'll
    # also need the first time step of the following day
    # Calculate the following day's date now, so we can be sure we have
    # the output directory
    t1 = datetime.datetime.strptime(yr+mth+day, "%Y%m%d")
    t2 = t1 + datetime.timedelta(days=1)
    yr2  = str(t2.year)
    mth2 = str(t2.month).zfill(2)
    day2 = str(t2.day).zfill(2)
    ecsdir='/'+yr+'/'+mth+'/'+day
    ecsdir2='/'+yr2+'/'+mth2+'/'+day2
    if not os.access(ecmwf_out+ecsdir, os.F_OK):
        os.makedirs(ecmwf_out+ecsdir)
    if not os.access(ecmwf_out+ecsdir2, os.F_OK):
        os.makedirs(ecmwf_out+ecsdir2)
    # If we have new ECMWF data available in the ECMWF NRT archive, do
    # the conversion into netcdf
    # Check what processed ECMWF files we already have
    proecm=os.listdir(ecmwf_out+'/'+yr+'/'+mth+'/'+day)
    # What data is available in the NRT archive
    rawecm=os.listdir(ecmwf_in_1+'/'+yr+'/'+mth+'/'+day)
    # We only want to keep the latest forecast file for each time slot. 
    # The production date/time is stored in characters 3-10 in the ECMWF
    # filename; characters 11-18 contain the target date-time, which are
    # at three hour intervals...
    tartimes=[mth+day+'000', mth+day+'030', mth+day+'060', 
              mth+day+'090', mth+day+'120', mth+day+'150', 
              mth+day+'180', mth+day+'210']
    jids = []
    ecm_outs = []
    for tartime in tartimes:
        if  abs(int(hr) - int(tartime[4:6])) <= 3:
            jid, ecm_out = process_latest_ecmwf_files(tartime, rawecm, proecm, ecmwf_in_1, 
                                 ecmwf_out, ecsdir, args)
            if jid and ecm_out:
                jids.append(jid)
                ecm_outs.append(ecm_out)
    # Now, if the processing is being run after the final ECMWF time 
    # slot for the day (9 pm), then we need to make sure we have the
    # first time slot from the following day available as well. So
    # repeat the above process for this timeslot           
    if int(hr) >= 21:
        proecm=os.listdir(ecmwf_out+'/'+yr2+'/'+mth2+'/'+day2)
        rawecm=os.listdir(ecmwf_in_1+'/'+yr2+'/'+mth2+'/'+day2)
        tartime=mth2+day2+'000'
        jid, ecm_path = process_latest_ecmwf_files(tartime, rawecm, proecm, ecmwf_in_1, 
                             ecmwf_out, ecsdir2, args)
        if jid and ecm_out:
                jids = np.append(jid)
                ecm_outs = np.append(ecm_out)
    return jids, ecm_outs


def process_pre(args, log_path, dependency=None, tag='pre'):
    """Call sequence for pre processor"""
    from pyorac.drivers import build_preproc_driver

    args = oracarg.check_args_preproc(args)
    driver, dependency = build_preproc_driver(args)

    # This must be called after building the driver as revision is unknown
    job_name = args.File.job_name(args.revision, tag)
    root_name = args.File.root_name(args.revision, args.processor, args.project,
                                    args.product_name)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, defaults.DIR_PERMISSIONS)

    out_file = os.path.join(args.out_dir, root_name + '.config.nc')
    if args.clobber >= CLOBBER['pre'] or not os.path.isfile(out_file):
        # Settings for batch processing
        values = {'job_name': job_name,
                  'log_file': os.path.join(log_path, job_name + '.log'),
                  'err_file': os.path.join(log_path, job_name + '.err')}
        if 'ram' in defaults.BATCH_VALUES:
                values['ram'] = defaults.BATCH_VALUES['ram'][1]
        else:
            values['ram'] = args.ram[1]
        if 'duration' in defaults.BATCH_VALUES:
            values['duration'] = defaults.BATCH_VALUES['duration'][1]
        else:
            values['duration'] = args.dur[1]
        if dependency is not None:
            values['depend'] = dependency
        exe = os.path.join(args.orac_dir, 'pre_processing', 'orac_preproc')
        if not os.path.isfile(exe):
            exe = os.path.join(args.orac_dir, 'orac_preproc')
        jid = call_exe(args, exe, driver, values)

    else:
        jid = None

    return jid, out_file


def process_main(args, log_path, tag='', dependency=None):
    """Call sequence for main processor"""
    from pyorac.drivers import build_main_driver

    args = oracarg.check_args_main(args)
    if args.phase == 'None':
        _, _, phase, _ = defaults.LUT_LOOKUP[args.lut_name](args.File, True)
        lutphs = args.lut_name
        if args.multilayer is not None:
            _, _, phase2, _ = defaults.LUT_LOOKUP[args.multilayer[0]](args.File, False)
            phase += "_" + phase2
    elif args.phase != 'CDF':
        lutphs = args.phase
        if args.multilayer is not None:
            phase = defin.SETTINGS[args.phase].name + "_" + defin.SETTINGS[args.multilayer[0]].name
        else:
            phase = defin.SETTINGS[args.phase].name

    job_name = args.File.job_name(tag=phase + tag)
    root_name = args.File.root_name(args.revision)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, defaults.DIR_PERMISSIONS)

    out_file = os.path.join(args.out_dir, root_name + lutphs + '.primary.nc')
    if args.clobber >= CLOBBER['main'] or not os.path.isfile(out_file):
        # Settings for batch processing
        values = {'job_name': job_name,
                  'log_file': os.path.join(log_path, job_name + '.log'),
                  'err_file': os.path.join(log_path, job_name + '.err')}
        if 'ram' in defaults.BATCH_VALUES:
                values['ram'] = defaults.BATCH_VALUES['ram'][2]
        else:
            values['ram'] = args.ram[2]
        if 'duration' in defaults.BATCH_VALUES:
            values['duration'] = defaults.BATCH_VALUES['duration'][2]
        else:
            values['duration'] = args.dur[2]
        if dependency is not None:
            values['depend'] = dependency

        driver = build_main_driver(args)
        exe = os.path.join(args.orac_dir, 'src', 'orac')
        if not os.path.isfile(exe):
            exe = os.path.join(args.orac_dir, 'orac')
        jid = call_exe(args, exe, driver, values)
    else:
        jid = None

    return jid, out_file


def process_post(args, log_path, files=None, dependency=None, tag='post'):
    """Call sequence for post processor"""
    from pyorac.drivers import build_postproc_driver

    args = oracarg.check_args_postproc(args)
    job_name = args.File.job_name(args.revision, tag)
    root_name = args.File.root_name(args.revision)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, defaults.DIR_PERMISSIONS)

    if files is None:
        # Find all primary files of requested phases in given input folders.
        files = []
        for phs in set(args.lut_names):
            for fdr in args.in_dir:
                files.extend(glob(os.path.join(
                    fdr, root_name + phs + '.primary.nc'
                )))

    if len(files) < 2:
        raise defin.FileMissing('sufficient processed files', args.target)

    out_file = os.path.join(
        args.out_dir, '.'.join(filter(
            None, (root_name, args.suffix, 'primary', 'nc')
        ))
    )
    if args.clobber >= CLOBBER['post'] or not os.path.isfile(out_file):
        # Settings for batch processing
        values = {'job_name': job_name,
                  'log_file': os.path.join(log_path, job_name + '.log'),
                  'err_file': os.path.join(log_path, job_name + '.err')}
        if 'ram' in defaults.BATCH_VALUES:
                values['ram'] = defaults.BATCH_VALUES['ram'][3]
        else:
            values['ram'] = args.ram[3]
        if 'duration' in defaults.BATCH_VALUES:
            values['duration'] = defaults.BATCH_VALUES['duration'][3]
        else:
            values['duration'] = args.dur[3]
        if dependency is not None:
            values['depend'] = dependency

        args.target = out_file
        driver = build_postproc_driver(args, files)
        exe = os.path.join(args.orac_dir, 'post_processing', 'orac_postproc')
        if not os.path.isfile(exe):
            exe = os.path.join(args.orac_dir, 'orac_postproc')
        jid = call_exe(args, exe, driver, values)

    else:
        jid = None

    return jid, out_file

def process_flux(args, log_path, files=None, dependency=None, tag='flux'):
    """Call sequence for post processor"""
    from glob import glob
    from pyorac.definitions import FileMissing, SETTINGS
    from pyorac.local_defaults import DIR_PERMISSIONS

    args = check_args_fluxes(args)
    job_name = args.File.job_name(args.revision, tag)
    root_name = args.File.root_name(args.revision)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, DIR_PERMISSIONS)

    if files is None:
        # Find all primary files of requested phases in given input folders.
        files = []
        for fdr in args.in_dir:
            out_dir_pre = args.out_dir + '/pre'
            files.extend(glob(os.path.join(
                    args.out_dir, root_name + '.primary.nc'
            )))
            files.extend(glob(os.path.join(
                    out_dir_pre, root_name + '.prtm.nc'
            )))
            files.extend(glob(os.path.join(
                    out_dir_pre, root_name + '.alb.nc'
            )))
            files.extend(glob(os.path.join(
                    out_dir_pre, root_name + '.config.nc'
            )))

    if len(files) < 4:
        raise FileMissing('sufficient processed files', args.target)
    out_file = os.path.join(
        args.out_dir, '.'.join(filter(
            None, (root_name, 'bugsrad', 'nc')
        ))
    )
    
    args.target = out_file
    
    if args.clobber >= CLOBBER['flux'] or not os.path.isfile(out_file):
        
        exe = os.path.join(args.orac_dir, '/derived_products/broadband_fluxes', 'process_broadband_fluxes')
        if not os.path.isfile(exe):
            exe = args.orac_dir+'/derived_products/broadband_fluxes/process_broadband_fluxes'
            
        
        cmd =exe + ' ' + files[0] + ' ' +files[1] + ' ' + files[2]+ ' ' + args.tsi+ ' ' + files[3]+ ' ' +args.target + \
                  ' ' + str(args.flux_alg)  + ' 0 0 0 0'
    
        if args.cci_aerpix:
            cmd += "cci_aerpix= " + files[0]
        
        if not args.batch:
            try:
                os.system(cmd)
                jid = None
            except CalledProcessError as err:
                raise OracError('{:s} failed with error code {:d}. {}'.format(
                    ' '.join(err.cmd), err.returncode, err.output
                ))
    
        else:
            # Write temporary script to call executable
            (gd, script_file) = tempfile.mkstemp('.sh', os.path.basename(exe)+'.',
                                             args.out_dir, True)
            g = os.fdopen(gd, "w")
            g.write("#!/bin/bash\n")
            # Define processing environment
            libs = read_orac_library_file(args.orac_lib)
            g.write("export LD_LIBRARY_PATH=" +
                              build_orac_library_path(libs) + "\n")
            g.write("export OPENBLAS_NUM_THREADS=1\n")
            try:
                g.write("export PPDIR=" + args.emos_dir + "\n")
            except AttributeError:
                pass
            BATCH.add_openmp_to_script(g)
            
            g.write(cmd+"\n")
            g.write("rm -f "+script_file+"\n")
            g.close()
            os.chmod(script_file, 0o700)
    
            try:
                # Collect batch settings from defaults, command line, and script
                batch_params = BATCH_VALUES.copy()
                batch_params['job_name'] = job_name
                batch_params['log_file'] = os.path.join(log_path, job_name + '.log')
                batch_params['err_file'] = os.path.join(log_path, job_name + '.err')
                batch_params['duration'] = '01:30:00'
                batch_params['ram'] = '3G'
                batch_params['procs'] = 1
                if values:
                    batch_params.update(values)
                batch_params.update({key: val for key, val in args.batch_settings})
    
                batch_params['procs'] = args.procs
    
                # Form batch queue command and call batch queuing system
                cmd = BATCH.list_batch(batch_params, exe=script_file)
    
                if args.verbose or args.script_verbose:
                    colour_print(' '.join(cmd), COLOURING['header'])
                out = check_output(cmd.split(' '), universal_newlines=True)
    
                # Parse job ID # and return it to the caller
                jid = BATCH.parse_out(out, 'ID')
            except CalledProcessError as err:
                raise OracError('Failed to queue job ' + exe)
            except SyntaxError as err:
                raise OracError(str(err))
            
    
    else:
        jid = None

    return jid, out_file

def process_flux(args, log_path, files=None, dependency=None, tag='flux'):
    """Call sequence for post processor"""

    args = oracarg.check_args_fluxes(args)
    job_name = args.File.job_name(args.revision, tag)
    root_name = args.File.root_name(args.revision)

    if not os.path.isdir(args.out_dir):
        os.makedirs(args.out_dir, defaults.DIR_PERMISSIONS)
    if files is None:
        # Find all primary files of requested phases in given input folders.
        files = []
        files.extend(glob(os.path.join(
                    args.out_dir, root_name + '.primary.nc'
        )))
        files.extend(glob(os.path.join(
                    args.pre_dir, root_name + '.prtm.nc'
        )))
        files.extend(glob(os.path.join(
                    args.pre_dir, root_name + '.alb.nc'
        )))
        files.extend(glob(os.path.join(
                    out_dir_pre, root_name + '.config.nc'
            )))

    if len(files) < 3:
        raise defin.FileMissing('sufficient processed files', args.target)
    out_file = os.path.join(
        args.out_dir, '.'.join(filter(
            None, (root_name, 'bugsrad', 'nc')
        ))
    )
    args.target = out_file
    
    if args.clobber >= CLOBBER['flux'] or not os.path.isfile(out_file):
        
        exe = os.path.join(args.orac_dir, '/derived_products/broadband_fluxes', 'process_broadband_fluxes')
        if not os.path.isfile(exe):
            exe = args.orac_dir+'/derived_products/broadband_fluxes/process_broadband_fluxes'
            
        
        cmd =exe + ' ' + files[0] + ' ' +files[1] + ' ' + files[2]+ ' ' + args.tsi+ ' ' + files[3]+ ' ' +args.target + \
                  ' ' + str(args.flux_alg)  + ' 0 0 0 0'
        if args.cci_aerpix:
            cmd += " cci_aerpix=" + files[0]
        if args.surface_to_process:
            cmd += " surface_to_process=" + str(args.surface_to_process)
        if args.verbose:
            cmd += " verbose=" + str(1)
        if args.procs:
            os.environ['OMP_NUM_THREADS'] = str(args.procs)
        if not args.batch:
            try:
                os.system(cmd)
                jid = None
            except CalledProcessError as err:
                raise OracError('{:s} failed with error code {:d}. {}'.format(
                    ' '.join(err.cmd), err.returncode, err.output
                ))
    
        else:
            # Write temporary script to call executable
            (gd, script_file) = tempfile.mkstemp('.sh', os.path.basename(exe)+'.',
                                             args.out_dir, True)
            g = os.fdopen(gd, "w")
            g.write("#!/bin/bash\n")
            # Define processing environment
            libs = read_orac_library_file(args.orac_lib)
            g.write("export LD_LIBRARY_PATH=" +
                              build_orac_library_path(libs) + "\n")
            g.write("export OPENBLAS_NUM_THREADS=1\n")
            try:
                g.write("export PPDIR=" + args.emos_dir + "\n")
            except AttributeError:
                pass
            defaults.BATCH.add_openmp_to_script(g)
            
            g.write(cmd+"\n")
            g.write("rm -f "+script_file+"\n")
            g.close()
            os.chmod(script_file, 0o700)
    
            try:
                # Collect batch settings from defaults, command line, and script
                batch_params = defaults.BATCH_VALUES.copy()
                batch_params['job_name'] = job_name
                batch_params['log_file'] = os.path.join(log_path, job_name + '.log')
                batch_params['err_file'] = os.path.join(log_path, job_name + '.err')
                batch_params['procs'] = 1
                if 'ram' in defaults.BATCH_VALUES:
                    batch_params['ram'] = defaults.BATCH_VALUES['ram'][4]
                else:
                    batch_params['ram'] = args.ram[4]
                if 'duration' in defaults.BATCH_VALUES:
                    batch_params['duration'] = defaults.BATCH_VALUES['duration'][4]
                else:
                    batch_params['duration'] = args.dur[4]
                if dependency is not None:
                    batch_params['depend'] = dependency
                batch_params.update({key: val for key, val in args.batch_settings})
    
                #batch_params['procs'] = args.procs
                # Form batch queue command and call batch queuing system
                cmd = defaults.BATCH.list_batch(batch_params, exe=script_file)
                out = check_output(cmd, universal_newlines=True)
    
                # Parse job ID # and return it to the caller
                jid = defaults.BATCH.parse_out(out, 'ID')
            except CalledProcessError as err:
                raise defin.OracError('Failed to queue job ' + exe)
            except SyntaxError as err:
                raise defin.OracError(str(err))
    
    else:
        jid = None

    return jid, out_file


def process_sisem_post(args, log_path, files=None, dependency=None, tag='sisem'):
    final_vars = ['time','lat','lon',
                  'boa_swdn_tot','boa_par_tot',
                  'boa_par_dif','boa_swdn_dif','stemp']
    job_name = args.File.job_name(args.revision, tag)
    root_name = args.File.root_name(args.revision)
    if files is None:
        # Find all primary files of requested phases in given input folders.
        files = glob(os.path.join(
                    args.out_dir, root_name + '.bugsrad.nc'))
    exe = args.orac_dir+'/tools/pyorac/sisem_postproc.py'
    cmd = 'python ' + exe + ' ' \
              + '--msi_root=' + root_name +' '+'--bugsradfile='+files[0]
    
    out_file = os.path.join(
        args.out_dir, '.'.join(filter(
            None, (root_name, 'sisem', 'nc'))))
    if not args.batch:
        try:
            os.system(cmd)
            jid = None
        except CalledProcessError as err:
            raise defin.OracError('{:s} failed with error code {:d}. {}'.format(
                ' '.join(err.cmd), err.returncode, err.output
            ))
    else:
        # Collect batch settings from defaults, command line, and script
        batch_params = defaults.BATCH_VALUES.copy()
        batch_params['job_name'] = job_name
        batch_params['log_file'] = os.path.join(log_path, job_name + '.log')
        batch_params['err_file'] = os.path.join(log_path, job_name + '.err')
        batch_params['procs'] = 1
        if 'ram' in defaults.BATCH_VALUES:
            batch_params['ram'] = defaults.BATCH_VALUES['ram'][5]
        else:
            batch_params['ram'] = args.ram[5]
        if 'duration' in defaults.BATCH_VALUES:
            batch_params['duration'] = defaults.BATCH_VALUES['duration'][5]
        else:
            batch_params['duration'] = args.dur[5]
        if dependency is not None:
                    batch_params['depend'] = dependency    
        batch_params.update({key: val for key, val in args.batch_settings})
        # Write temporary script to call executable
        
        (gd, script_file) = tempfile.mkstemp('.sh', os.path.basename(exe)+'.',
                                             args.out_dir, True)
        g = os.fdopen(gd, "w")
        g.write("#!/bin/bash\n")
        # Define processing environment
        libs = read_orac_library_file(args.orac_lib)
        g.write("source /home/users/$USER/miniforge3/bin/activate\n")
        g.write("conda activate sev_ml_tf\n")
        g.write("export LD_LIBRARY_PATH=" +
                              build_orac_library_path(libs) + "\n")
        g.write("export OPENBLAS_NUM_THREADS=1\n")
        defaults.BATCH.add_openmp_to_script(g) 
        g.write(cmd+"\n")
        g.write("rm -f "+script_file+"\n")
        g.close()
        os.chmod(script_file, 0o700)
        cmd = defaults.BATCH.list_batch(batch_params, exe=script_file)
        if args.verbose or args.script_verbose:
                    colour_print(' '.join(cmd), defin.COLOURING['header'])
        out = check_output(cmd, universal_newlines=True)
        # Parse job ID # and return it to the caller
        jid = defaults.BATCH.parse_out(out, 'ID')
    return jid, out_file


def call_reformat(args, log_path, exe, out_file, dependency=None):
    """Reformat outputs using the script provided."""
    from pyorac.colour_print import colour_print
    from subprocess import check_call, check_output, CalledProcessError

    # Optionally print command and driver file contents to StdOut
    if args.verbose or args.script_verbose or args.dry_run:
        colour_print('{} {} 1 <<<'.format(exe, out_file), defin.COLOURING['header'])

    if args.dry_run:
        return -1

    job_name = args.File.job_name(args.revision, 'format')

    if not args.batch:
        try:
            check_call([exe, out_file, "1"])
        except CalledProcessError as err:
            raise defin.OracError('{:s} failed with error code {:d}. {}'.format(
                ' '.join(err.cmd), err.returncode, err.output
            ))

    else:
        try:
            # Collect batch settings from defaults, command line, and script
            batch_params = defaults.BATCH_VALUES.copy()
            batch_params['job_name'] = job_name
            batch_params['log_file'] = os.path.join(log_path, job_name + '.log')
            batch_params['err_file'] = os.path.join(log_path, job_name + '.err')
            batch_params['duration'] = '01:00'
            batch_params['ram'] = 5000
            batch_params['procs'] = 1
            if dependency is not None:
                batch_params['depend'] = dependency
            batch_params.update({key: val for key, val in args.batch_settings})

            # Form batch queue command and call batch queuing system
            cmd = defaults.BATCH.list_batch(batch_params, exe=[exe, out_file, "1"])

            if args.verbose or args.script_verbose:
                colour_print(' '.join(cmd), defin.COLOURING['header'])
            out = check_output(cmd, universal_newlines=True)

            # Parse job ID # and return it to the caller
            jid = defaults.BATCH.parse_out(out, 'ID')
            return jid
        except CalledProcessError:
            raise defin.OracError('Failed to queue job ' + exe)
        except SyntaxError as err:
            raise defin.OracError(str(err))


def process_all(orig_args):
    """Run the ORAC pre, main, and post processors on a file."""
    from argparse import ArgumentParser

    # Generate main-processor-only parser
    pars = ArgumentParser()
    oracarg.args_common(pars)
    oracarg.args_main(pars)
    # We need one argument from args_cc4cl()
    pars.add_argument("--sub_dir", default="")
    compare = pars.parse_args("")

    # Copy input arguments as we'll need to fiddle with them
    orig_args = oracarg.check_args_common(orig_args)
    orig_args = oracarg.check_args_cc4cl(orig_args)
    log_path = os.path.join(orig_args.out_dir, defaults.LOG_DIR)
    args = deepcopy(orig_args)

    written_dirs = set()  # The folders we actually wrote to

    # Work out output filename
    args.out_dir = os.path.join(orig_args.out_dir, defaults.PRE_DIR)

    jid_pre, _ = process_pre(args, log_path, tag="pre{}".format(args.label))
    if jid_pre is not None:
        written_dirs.add(args.out_dir)

    # Run main processor -------------------------------------------------------
    root_name = args.File.root_name(args.revision, args.processor, args.project,
                                    args.product_name)
    args.target = root_name + ".config.nc"
    out_files = []  # All files that would be made (facilitates --dry_run)
    jid_main = []  # ID no. for each queued job
    args.in_dir = [args.out_dir]
    jid = None
    for sett in args.settings:
        phs_args = deepcopy(args)
        parsed_settings_arguments = pars.parse_args(sett.split())
        for key, val in compare.__dict__.items():
            if val == parsed_settings_arguments.__dict__[key]:
                parsed_settings_arguments.__dict__.pop(key)
        phs_args.__dict__.update(parsed_settings_arguments.__dict__)
        phs_args.out_dir = os.path.join(orig_args.out_dir, phs_args.sub_dir)

        jid, out = process_main(phs_args, log_path, dependency=jid_pre,
                                tag=phs_args.sub_dir + phs_args.label)
        out_files.append(out)
        if jid is not None:
            jid_main.append(jid)
            written_dirs.add(args.out_dir)

    # Run postprocessor if necessary
    if len(args.settings) > 1:
        args.target = root_name + "NULL.primary.nc"
        args.in_dir = written_dirs
        args.out_dir = orig_args.out_dir
        jid, out_file = process_post(
            args, log_path, out_files, dependency=jid_main,
            tag="post{}".format(args.label)
        )
        if jid is not None:
            written_dirs.add(args.out_dir)
    else:
        out_file = out_files[0]

    # Run CCI formatting
    if args.reformat != "":
        call_reformat(args, log_path, args.reformat, out_file, dependency=jid)

    # Output root filename and output folders for regression tests
    return jid, out_file


def run_regression(in_file):
    """Run the regression test on a set of ORAC files."""
    import re

    from pyorac.util import compare_orac_out
    from warnings import warn

    regex = re.compile(r"_R(\d+)")
    this_revision = int(in_file.revision)
    for fdr in in_file.folders:
        for this_file in glob(os.path.join(
                fdr, "**", in_file.root_name() + "*nc"
        ), recursive=True):
            # Find previous file version
            old_revision = 0
            old_file = None
            for filename in glob(regex.sub("_R[0-9][0-9][0-9][0-9]", this_file)):
                rev = int(regex.search(filename).group(1))
                if old_revision < rev < this_revision:
                    old_revision = deepcopy(rev)
                    old_file = deepcopy(filename)

            if old_file is None:
                warn("Could not locate previous file: " + this_file, defin.OracWarning)
                continue

            compare_orac_out(this_file, old_file)
