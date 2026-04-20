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
#   run radialplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0 --multi_panel
#   run radialplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
#   run radialplot.py --mode lsmzsfr --lsm 9.5,10.5
#   run radialplot.py --quant electron --mode lsmzsfr --lsm all
#   run radialplot.py --quant electron --mode lsmzsfr --lsm all --multi_panel
#   run radialplot.py --quant gas --mode lsmzsfr --lsm all
#   run radialplot.py --quant gas --mode lsmzsfr --lsm all --multi_panel

#	--------------------------	Import modules	---------------------------
from craft_utils import *
setup_plot_style()
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

    if args.multi_panel:
        nrows, ncols = get_grid_size(len(df_snap))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
        fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

    for i, snap in df_snap.iterrows():
        thisfile = args.radial_data_dir / f'{snap["snap"]}_{snap["halo"]}_FRB{args.quant_text}_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_radprof.pkl'
        with open(thisfile, 'rb') as file_obj: cube	= pkl.load(file_obj)

        title	= r"log ($M_* / M_{\odot}$) = %.2f, SFR = %.2f $M_{\odot} yr^{-1}$"%(snap["lsm"], snap["sfr"])
        outfile = args.radial_plot_dir / f'{args.mode}_{snap["halo"]}_{snap["snap"]}_{args.quant}_density_radprof.pdf'
        ax = pfns.plot_nerad([cube], incranges, title, outfile, 3.0, hide=args.hide, given_ax=axes[i // ncols][i % ncols] if args.multi_panel else None, fortalk=args.fortalk)

        if args.multi_panel:
            if i // ncols < nrows - 1:
                ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax.set_xlabel('')
            if i % ncols > 0:
                ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax.set_ylabel('')
    
    if args.multi_panel: save_fig(fig, args.fig_dir, f'{args.mode}_{args.quant}_density_multipanel_radprof.pdf', args)
    return

# -----------------------------------------------------------------------------
def execute_mode_halo(df_snap, incranges, args):
    '''
    Function to execute mode halo
    Returns nothing
    '''
    print(f"\nTracking halo {args.halo} ...\n")
    df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]
    
    if args.multi_panel:
        nrows, ncols = get_grid_size(len(df_snap))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
        fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

    for i, snap in df_snap.iterrows():
        thisfile = args.radial_data_dir / f'{snap["snap"]}_{snap["halo"]}_FRB{args.quant_text}_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_radprof.pkl'
        with open(thisfile, 'rb') as file_obj: cube	= pkl.load(file_obj)

        title	= r"log ($M_* / M_{\odot}$) = %.2f, SFR = %.2f $M_{\odot} yr^{-1}$"%(snap["lsm"], snap["sfr"])
        outfile = args.radial_plot_dir / f'{args.mode}_{snap["halo"]}_{snap["snap"]}_{args.quant}_density_radprof.pdf'
        ax = pfns.plot_nerad([cube], incranges, title, outfile, 3.0, hide=args.hide, given_ax=axes[i // ncols][i % ncols] if args.multi_panel else None, fortalk=args.fortalk)

        if args.multi_panel:
            if i // ncols < nrows - 1:
                ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                ax.set_xlabel('')
            if i % ncols > 0:
                ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                ax.set_ylabel('')

    if args.multi_panel: save_fig(fig, args.fig_dir, f'{args.mode}_{args.halo}_{args.quant}_density_multipanel_radprof.pdf', args)
    return

# -----------------------------------------------------------------------------
def execute_mode_lsmzsfr(df_snap, incranges, args, given_ax=None):
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

    ax = pfns.plot_nerad(cubes, incranges, title, outfile, 2.8, hide=args.hide, subtitle=subtitle, given_ax=given_ax, fortalk=args.fortalk)

    return ax

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    args = parse_args()
    if not args.keep: plt.close('all')
    args.radial_data_dir = args.data_dir / f'radial_profiles_{args.quant}_density'
    args.radial_plot_dir = args.plot_dir / f'{args.quant}_density_profiles'
    args.quant_text = '_El_number' if args.quant == 'electron' else ''
    incranges	=	np.array([[item - dinc/2, item + dinc/2] for item in incvals])

    # ------------setup multi-panel figure if needed-------
    if args.multi_panel and args.mode == 'lsmzsfr':
        nrows, ncols = get_grid_size(len(args.lsm_bins) * len(args.lsfr_bins))
        fig, axes = plt.subplots(nrows, ncols, figsize=(10, 8))
        fig.subplots_adjust(left=0.07, bottom=0.07, right=0.98, top=0.98, wspace=0.01, hspace=0.01)

    # ------------looping over SFR bins----------
    for index2, this_lsfr_bin in enumerate(args.lsfr_bins):
        print(f'\n\tRunning ({index2 + 1}/{len(args.lsfr_bins)}) for SFR bin {this_lsfr_bin}..\n')
        args.lsfr_range = this_lsfr_bin

        # ------------looping over stellar mass bins----------
        for index, this_lsm_bin in enumerate(args.lsm_bins):
            print(f'\n\t\tRunning ({index + 1}/{len(args.lsm_bins)}) for stellar mass bin {this_lsm_bin}..\n')
            args.lsm_range = this_lsm_bin

            if args.multi_panel:
                nrow = index if len(args.lsfr_bins) > 1 else index // ncols
                ncol = index2 if len(args.lsfr_bins) > 1 else index % ncols

            #	-------------------------	Initialize	-----------------------------------
            filespecs =	args.data_dir / "lsm_sfr.txt"
            df_snap = pd.read_csv(filespecs, sep=r'\s+', engine='python')
            df_snap = df_snap.rename(columns={f'{df_snap.columns[0]}':f'{df_snap.columns[0][1:]}'})
            df_snap['log_sfr'] = np.log10(df_snap['sfr'])
            df_snap = df_snap[(df_snap['redshift'].between(args.z_range[0], args.z_range[1])) & 
                            (df_snap['log_star_mass'].between(args.lsm_range[0], args.lsm_range[1])) & 
                            (df_snap['log_sfr'].between(args.lsfr_range[0], args.lsfr_range[1]))].reset_index(drop=True)
            print (f'\t\tFound {len(df_snap)} snapshots, within log mass range {args.lsm_range}, log sfr range {args.lsfr_range} and redshift range {args.z_range}')
            if (len(df_snap) < 1): sys.exit('Exiting because no snapshot found... ')
            
            #	-------------------------	Execute tasks	-------------------------------
            if (args.mode=='indi'):
                execute_mode_indi(df_snap, incranges, args)
            elif (args.mode=='halo'):
                execute_mode_halo(df_snap, incranges, args)            
            elif (args.mode=='lsmzsfr'):
                # ---------------make the plots-----------------
                ax = execute_mode_lsmzsfr(df_snap, incranges, args, given_ax=axes[nrow][ncol] if args.multi_panel else None)

                if args.multi_panel:
                    if nrow < nrows - 1:
                        ax.tick_params(axis='x', which='major', labelsize=0, labelbottom=False)
                        ax.set_xlabel('')
                    if ncol > 0:
                        ax.tick_params(axis='y', which='major', labelsize=0, labelbottom=False)
                        ax.set_ylabel('')
            else:
                print("\n\t\tHmm...What mode is that again...?\n")

    if args.mode == 'lsmzsfr' and args.multi_panel:
        save_fig(fig, args.fig_dir, f'{args.mode}_inc_{args.inc_range[0]}_{args.inc_range[1]}_{args.quant}_density_multipanel_radprof.pdf', args)

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































