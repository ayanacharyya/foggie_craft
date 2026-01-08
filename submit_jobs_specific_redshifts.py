##!/usr/bin/env python3

"""

    Title :      submit_jobs
    Notes :      Python wrapper to create and submit one or more jobs on pleiades
    Author :     Ayan Acharyya
    Started :    July 2021
    Example :    run submit_jobs_specific_redshifts.py --call make_3D_FRB_electron_density --system ayan_pleiades --do_all_halos --queue ldan --mem 1500GB --prefix cmzs --do_redshifts 4,3,2,1.5,1,0.8,0.6,0.5,0.4,0.3,0.2,0.1,0 --opt_args "--res 1 --upto_kpc 10 --docomoving --use_cen_smoothed"
"""
import subprocess, argparse, datetime, os
from collections import defaultdict
import numpy as np
import pandas as pd

# ------------------------------------------------------
def execute_command(command, is_dry_run=False):
    '''
    Function to decide whether to execute a command or simply print it out (for dry run)
    '''
    if is_dry_run:
        print('Not executing command:', command, '\n')
        return -99
    else:
        print('Executing command:', command, '\n')
        job = subprocess.check_output([command], shell=True)[:-1]
        return job

# ---------------------------------------------------------
def parse_args():
    '''
    Function to parse keyword arguments
    '''
    parser = argparse.ArgumentParser(description="calling plotobservables for full parameter space")
    parser.add_argument('--system', metavar='system', type=str, action='store', default='ayan_pleiades')
    parser.add_argument('--proj', metavar='proj', type=str, action='store', default='s2358')
    parser.add_argument('--queue', metavar='queue', type=str, action='store', default='long')
    parser.add_argument('--nnodes', metavar='nnodes', type=int, action='store', default=1)
    parser.add_argument('--ncores', metavar='ncores', type=int, action='store', default=None)
    parser.add_argument('--nhours', metavar='nhours', type=int, action='store', default=None)
    parser.add_argument('--ncpus', metavar='ncpus', type=int, action='store', default=None)
    parser.add_argument('--nmins', metavar='nmins', type=int, action='store', default=0)
    parser.add_argument('--proc', metavar='proc', type=str, action='store', default='has')
    parser.add_argument('--memory', metavar='memory', type=str, action='store', default=None)
    parser.add_argument('--aoe', metavar='aoe', type=str, action='store', default=None)
    parser.add_argument('--start', metavar='start', type=int, action='store', default=1)
    parser.add_argument('--stop', metavar='stop', type=int, action='store', default=1)
    parser.add_argument('--halo', metavar='halo', type=str, action='store', default=None)
    parser.add_argument('--run', metavar='run', type=str, action='store', default='nref11c_nref9f')
    parser.add_argument('--prefix', metavar='prefix', type=str, action='store', default=None)
    parser.add_argument('--callfunc', metavar='callfunc', type=str, action='store', default='filter_star_properties')
    parser.add_argument('--dryrun', dest='dryrun', action='store_true', default=False)
    parser.add_argument('--do_all_sims', dest='do_all_sims', action='store_true', default=False)
    parser.add_argument('--do_all_halos', dest='do_all_halos', action='store_true', default=False)
    parser.add_argument('--galrad', metavar='galrad', type=str, action='store', default=None)
    parser.add_argument('--snapstart', metavar='snapstart', type=int, action='store', default=30)
    parser.add_argument('--snapstop', metavar='snapstop', type=int, action='store', default=30)
    parser.add_argument('--do_redshifts', metavar='do_redshifts', type=str, action='store', default=None)
    parser.add_argument('--opt_args', metavar='opt_args', type=str, action='store', default='')
    args, leftovers = parser.parse_known_args()

    return args

# ---------------------------------------------------------
if __name__ == '__main__':
    time_of_begin = datetime.datetime.now()
    args = parse_args()
    if args.do_redshifts is not None: args.do_redshifts = [float(item) for item in args.do_redshifts.split(',')]
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'

    # ----------special settings for ldan queue--------
    if args.queue == 'ldan':
        args.proc = 'ldan'
        args.nnodes = 1
        args.ncores = None
    # ----------special settings for endeavour queue--------
    elif args.queue[:2] == 'e_':
        args.proc = 'cas_end'

    #----------setting different variables based on args--------
    systemflag = ' --system ' + args.system
    runsimflag = ' --run ' + args.run
    prefixtext = 'frb_'

    if 'pleiades' in args.system: jobscript_path = '/nobackupp19/aachary2/foggie_craft/foggie_craft/'
    elif args.system == 'ayan_local': jobscript_path = os.getenv('HOME') + '/Work/astro/ayan_codes/foggie_craft/'

    if args.system == 'ayan_local': jobscript_template = 'jobscript_template_ayan_pleiades.txt'
    else: jobscript_template = 'jobscript_template_' + args.system + '.txt'

    callfile = jobscript_path + args.callfunc + '.py'

    if 'pleiades' in args.system or args.system == 'ayan_local':
        procs_dir = {'san':(16, 32), 'ivy':(20, 64), 'has':(24, 128), 'bro':(28, 128), 'bro_ele':(28, 128), 'sky_ele':(40, 192), 'cas_ait':(40, 192), 'ldan':(16, 750), 'cas_end':(28, 185)} # (nnodes, mem) for each proc, from https://www.nas.nasa.gov/hecc/support/kb/pbs-resource-request-examples_188.html
        max_hours_dict = defaultdict(lambda: 120, low=4, normal=8, long=120, e_long=72, e_normal=8, e_vlong=600, e_debug=2, debug=2, devel=2, ldan=72) # from https://www.nas.nasa.gov/hecc/support/kb/pbs-job-queue-structure_187.html
        if 'pleiades' in args.system: workdir = '/nobackupp19/aachary2/foggie_craft/pleiades_workdir' # for pleiades
        elif args.system == 'ayan_local': workdir = '.'
        nnodes = args.nnodes
        ncores = args.ncores if args.ncores is not None else procs_dir[args.proc][0]
        memory = args.memory if args.memory is not None else str(procs_dir[args.proc][1]) + 'GB' # minimum memory per node; by default the entire node me is allocated, therefore it is redundant to specify mem as the highest available memory per node
        qname = args.queue
        if args.queue[:2] == 'e_':  qname += '@pbspl4' # add this for endeavour nodes

    os.chdir(workdir)

    # ----------determining what resource request goes into the job script, based on queues, procs, etc.---------
    nhours = args.nhours if args.nhours is not None else '01' if args.dryrun or args.queue == 'devel' else '%02d' % (max_hours_dict[args.queue])
    if args.ncpus is None:
        if args.do_redshifts is None: ncpus = nnodes * ncores
        else: ncpus = min(len(args.do_redshifts), nnodes * ncores)
    else:
        ncpus = args.ncpus

    resources = 'select=' + str(nnodes) + ':ncpus=' + str(ncores)

    if args.queue[:2] == 'e_': resources += ':mem=' + memory # for submitting to endeavour
    else: resources += ':mpiprocs=' + str(ncores)

    if args.queue == 'ldan': resources += ':mem=' + memory # may specify mem per node for jobs on LDAN (for other procs it is by default the max available node mem)
    else: resources += ':model=' + args.proc # need to specify the proc (but not necessarily the mem if I'm using the full node memory) if not an LDAN job

    if args.aoe is not None: resources += ':aoe=' + args.aoe

    #----------looping over and creating + submitting job files--------
    halos = ['8508', '5036', '5016', '4123', '2392', '2878']
    if args.do_all_halos: args.stop = len(halos) # do all halos by submitting ..

    for jobid in range(args.start, args.stop+1):
        thishalo = halos[jobid - 1] if args.halo is None else args.halo
        haloflag = ' --halo ' + thishalo
        jobname = prefixtext + thishalo
        if jobname[:3] != args.proc: jobname = args.proc + '_' + jobname

        # ---------determining which RD/DD outputs to run on----------------
        if args.do_redshifts is None:
            outputs = args.output
        else:
            df = pd.read_csv(args.code_dir + f'halo_infos/00{thishalo}/nref11c_nref9f/halo_cen_smoothed', sep=r'\s*\|\s*', engine='python')
            df = df.dropna(axis=1, how='all')[['snap', 'redshift']]
            output_list = []
            for redshift in args.do_redshifts:
                idx = (df['redshift'] - redshift).abs().idxmin()
                output_list.append(df.loc[idx, 'snap'])
            outputs = ','.join(output_list)
        print(f'Halo {thishalo}: Going to submit jobs for {len(outputs.split(","))} outputs: {outputs}..')

        # ----------replacing keywords in jobscript template to make the actual jobscript---------
        out_jobscript = workdir + '/jobscript_' + jobname + '.sh'

        replacements = {'PROJ_CODE': args.proj, 'RUN_NAME': jobname, 'NHOURS': nhours, 'NMINS': args.nmins, 'CALLFILE': callfile, 'WORKDIR': workdir, \
                        'QNAME': qname, 'RESOURCES': resources, 'RUNSIMFLAG': runsimflag, 'OUTPUT_FLAG': outputs, \
                        'SYSTEMFLAG': systemflag, 'NCPUS': str(ncpus), 'HALOFLAG': haloflag, 'NSECONDS':str(int(nhours) * 3600), 'OPT_ARGS': args.opt_args} # keywords to be replaced in template jobscript

        with open(jobscript_path + jobscript_template) as infile, open(out_jobscript, 'w') as outfile:
            for line in infile:
                for src, target in replacements.items():
                    line = line.replace(str(src), str(target))
                outfile.write(line) # replacing and creating new jobscript file

        print('Going to submit job ' + jobname+' at {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
        jobid = execute_command('qsub ' + out_jobscript, is_dry_run=args.dryrun)

    print('Submitted all '+str(args.stop - args.start + 1)+' job/s from index '+str(args.start)+' to '+str(args.stop) + '\n')