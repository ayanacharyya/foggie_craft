#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								AB, August 2024
#   Modified by AA in Apr 2026
#   
#   Examples of how to run (from within ipython):   
#   run radialplot.py --mode lsmzsfr --rangekpc 200 --reskpc 0.5 --z_range 0,6 --inc 0,90 --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run radialplot.py --mode lsmzsfr --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run radialplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0
#   run radialplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
#   run radialplot.py --mode lsmzsfr --lsm 9.5,10.5
#   run radialplot.py --quant electron --mode lsmzsfr --lsm all
#   run radialplot.py --quant gas --mode lsmzsfr --lsm all

#	--------------------------	Import modules	---------------------------
from craft_utils import *
import plotfns as pfns

start_time = datetime.now()

# -----------------------------------------------------------------------------
def print_instructions ():

	#	Print instructions to terminal
	
	print("\n            You probably need some assistance here!\n")
	print("\n Arguments are     --- <mode> <zmin/zmax> <args.lsm_range[0]/args.lsm_range[1]> <args.lsfr_range[0]/args.lsfr_range[1]> <rangkpc> <reskpc> <args.inc_range[0]/args.inc_range[1]> <halo> <figname> <args.resfile_prefix>\n")
	print(" Supported Modes are --- indi     (Plot individual files)")
	print("                     --- halo     (Track a particular halo)")
	print("                     --- lsmzsfr  (Combine all wihin the lsm and z range)")
	
	
	print("\n            Now let's try again!\n")
	
	return(0)

# -----------------------------------------------------------------------------
def execute_mode_indi(df_snap, incranges, args):
    '''
    Function to execute mode indi
    Returns nothing
    '''
    print("\nPloting individual snaps...\n")

    for i, snap in df_snap.iterrows():
        thisfile = args.radial_data_dir / f'{snap["snap"]}_{snap["halo"]}_FRB{args.quant_text}_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_radprof.pkl'
        with open(thisfile, 'rb') as file_obj: cube	= pkl.load(file_obj)

        title	= r"log ($M_* / M_{\odot}$) = %.2f, SFR = %.2f $M_{\odot} yr^{-1}$"%(snap["lsm"], snap["sfr"])
        outfile = args.radial_plot_dir / f'{args.mode}_{snap["halo"]}_{snap["snap"]}_radprof.pdf'
        pfns.plot_nerad([cube], incranges, title, outfile, 3.0, hide=args.hide)

    return

# -----------------------------------------------------------------------------
def execute_mode_halo(df_snap, incranges, args):
    '''
    Function to execute mode halo
    Returns nothing
    '''
    print(f"\nTracking halo {args.halo} ...\n")
    df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]
    
    for i, snap in df_snap.iterrows():
        thisfile = args.radial_data_dir / f'{snap["snap"]}_{snap["halo"]}_FRB{args.quant_text}_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_radprof.pkl'
        with open(thisfile, 'rb') as file_obj: cube	= pkl.load(file_obj)

        title	= r"log ($M_* / M_{\odot}$) = %.2f, SFR = %.2f $M_{\odot} yr^{-1}$"%(snap["lsm"], snap["sfr"])
        outfile = args.radial_plot_dir / f'{args.mode}_{snap["halo"]}_{snap["snap"]}_radprof.pdf'
        pfns.plot_nerad([cube], incranges, title, outfile, 3.0, hide=args.hide)

    return

# -----------------------------------------------------------------------------
def execute_mode_lsmzsfr(df_snap, incranges, args):
    '''
    Function to execute mode lsmzsfr
    Returns nothing
    '''
    print(f"\nCombining all within lsm range {args.lsm_range} and z range {args.z_range} \n")

    cubes = []

    for i, snap in df_snap.iterrows():
        thisfile = args.radial_data_dir / f'{snap["snap"]}_{snap["halo"]}_FRB{args.quant_text}_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_radprof.pkl'
        with open(thisfile, 'rb') as file_obj: cube	= pkl.load(file_obj)
        cubes.append(cube)

    print(f"Total number of cubes = {len(cubes)}")	

    median_lsm	= np.nanmedian(df_snap["log_star_mass"]) 
    median_sfr	= np.nanmedian(df_snap["sfr"])

    #title	= "log ($M_* / M_{\odot}$) = %.2f, SFR = %.2f $M_{\odot} yr^{-1}$"%(median_lsm, median_sfr)
    title	= "%.1f < log ($M_* / M_{\odot}$) < %.1f"%(args.lsm_range[0], args.lsm_range[1])
    subtitle	= "%.1f < log (SFR / $M_{\odot} yr^{-1}$) < %.1f"%(args.lsfr_range[0], args.lsfr_range[1])
    outfile = args.radial_plot_dir / f'{args.mode}_{snap["halo"]}_{snap["snap"]}_radprof.pdf'

    pfns.plot_nerad(cubes, incranges, title, outfile, 2.8, hide=args.hide, subtitle=subtitle)

    return

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    args = parse_args()
    args.radial_data_dir = args.data_dir / f'radial_profiles_{args.quant}_density'
    args.radial_plot_dir = args.plot_dir / f'{args.quant}_density_profiles'
    args.quant_text = '_El_number' if args.quant == 'electron' else ''
    incranges	=	np.array([[item - dinc/2, item + dinc/2] for item in incvals])

    # ------------looping over stellar mass bins----------
    for index, this_lsm_bin in enumerate(args.lsm_bins):
        print(f'\nRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
        args.lsm_range = this_lsm_bin
        #	-------------------------	Initialize	-----------------------------------
        filespecs =	args.data_dir / "lsm_sfr.txt"
        df_snap = pd.read_csv(filespecs, sep=r'\s+', engine='python')
        df_snap = df_snap.rename(columns={f'{df_snap.columns[0]}':f'{df_snap.columns[0][1:]}'})
        df_snap['log_sfr'] = np.log10(df_snap['sfr'])
        df_snap = df_snap[(df_snap['redshift'].between(args.z_range[0], args.z_range[1])) & 
                        (df_snap['log_star_mass'].between(args.lsm_range[0], args.lsm_range[1])) & 
                        (df_snap['log_sfr'].between(args.lsfr_range[0], args.lsfr_range[1]))]
        print (f'Found {len(df_snap)} snapshots, within log mass range {args.lsm_range}, log sfr range {args.lsfr_range} and redshift range {args.z_range}')
        if (len(df_snap) < 1): sys.exit('Exiting because no snapshot found... ')
        
        #	-------------------------	Execute tasks	-------------------------------
        if (args.mode=='indi'):
            execute_mode_indi(df_snap, incranges, args)
        elif (args.mode=='halo'):
            execute_mode_halo(df_snap, incranges, args)            
        elif (args.mode=='lsmzsfr'):
            execute_mode_lsmzsfr(df_snap, incranges, args)
        else:
            print("\nHmm...What mode is that again...?\n")

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































