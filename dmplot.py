#
#	Script for estimating FRB host galaxy DM from simulated electron density cubes
#
#								AB, August 2024
#   Modified by AA in Apr 2026
#   
#   Examples of how to run (from within ipython):   
#   run dmplot.py --mode lsmzsfr --rangekpc 200 --reskpc 0.5 --z_range 0,6 --inc 0,90 --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run dmplot.py --mode lsmzsfr --lsm 9.5,10.5 --lsfr 0,0.5 --resfile_prefix all_lsm
#   run dmplot.py --mode indi --lsm 9.5,10.0 --lsfr=-1,0
#   run dmplot.py --mode halo --halo 5036 --lsm 9.5,10.0 --lsfr=-1,0
#   run dmplot.py --mode lsmzsfr --lsm 9.5,10.5
#   run dmplot.py --mode lsmzsfr --lsm all

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
def execute_mode_indi(df_snap, args):
    '''
    Function to execute mode indi
    Returns nothing
    '''
    print("\nPloting individual snaps...\n")

    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        pars, epars	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 3.0, hide=args.hide, bin_col='impf', data_col='losdm')

        lgsm	= np.log10(10.0 ** snap["log_gas_mass"] + 10.0 ** snap["log_star_mass"])

        # --------------initialise dataframe------------------
        df_out = pd.DataFrame({'lsm_bin': pd.Interval(args.lsm_range[0], args.lsm_range[1]),
                            'lsfr_bin': pd.Interval(args.lsfr_range[0], args.lsfr_range[1]),
                            'medlsm': snap["log_star_mass"],
                            'medsfr': snap["sfr"],
                            'medlgm': snap["log_gas_mass"],
                            'medlgsm': lgsm,
                            'ledssfr9': np.log10(snap["sfr"]) - snap["log_star_mass"] + 9,
                            'medgsfr9': np.log10(snap["sfr"]) - snap["log_gas_mass"] + 9,
                            'medgssfr9': np.log10(snap["sfr"]) - lgsm + 9,
                            'r0': pars[0],
                            'er0': epars[0],
                            'D0': pars[1],
                            'eD0': epars[1],
                            }, index=[0])
        
        outfile = f'{args.resfile_prefix}_indiv_inc_{args.inc_range[0]}_{args.inc_range[1]}.txt'
        df_out.to_csv(outfile, mode='a', sep='\t', header=not os.path.exists(outfile), index=None)

    return

# -----------------------------------------------------------------------------
def execute_mode_halo(df_snap, args):
    '''
    Function to execute mode halo
    Returns nothing
    '''
    print(f"\nTracking halo {args.halo} ...\n")
    df_snap = df_snap[df_snap["halo"].astype(str) == args.halo]
    
    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]

        print(f'{snap["snap"]}_{snap["halo"]}: Total number of LoS = {len(this_df)}')
        print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

        outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{snap["halo"]}_{snap["snap"]}'
        
        pfns.pltdm_ind_imf(this_df, snap['log_star_mass'], snap['sfr'], args.inc_range, snap['redshift'], outfile, 3.0, hide=False, bin_col1='impf', bin_col2='distmaj', data_col='losdm')
        
        pars, epars	= pfns.pltdm_ind_imf_1d(this_df, snap['log_star_mass'], snap['sfr'], args.lsm_range, outfile + '_1d', 3.0, hide=args.hide, bin_col='impf', data_col='losdm')

    return

# -----------------------------------------------------------------------------
def execute_mode_lsmzsfr(df_snap, args):
    '''
    Function to execute mode lsmzsfr
    Returns nothing
    '''
    print(f"\nCombining all within lsm range {args.lsm_range} and z range {args.z_range} \n")

    combined_df = pd.DataFrame()
    
    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        dm_arr	= np.load(thisfile)
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df = this_df[this_df['inc'].between(args.inc_range[0], args.inc_range[1])]
        combined_df = pd.concat([combined_df, this_df], ignore_index=True)
        
    print(f"Total number of LoS = {len(combined_df)}")	
    print("Plotting DMs within inclination ",args.inc_range[0], args.inc_range[1])

    median_lsm = np.nanmedian(df_snap["log_star_mass"])
    median_lgm = np.nanmedian(df_snap["log_gas_mass"])
    median_sfr = np.nanmedian(df_snap["sfr"])
    lgsm	= np.log10(10.0 ** df_snap["log_gas_mass"] + 10.0 ** df_snap["log_star_mass"])
    
    outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}/{args.mode}_lsm_{args.lsm_range[0]}_{args.lsm_range[1]}_lsfr_{args.lsfr_range[0]}_{args.lsfr_range[1]}_1d'
    pars, epars	= pfns.pltdm_ind_imf_1d(combined_df, median_lsm, median_sfr, args.lsm_range, outfile, 2.6, hide=args.hide, bin_col='impf', data_col='losdm')

    # --------------initialise dataframe------------------
    df_out = pd.DataFrame({'lsm_bin': pd.Interval(args.lsm_range[0], args.lsm_range[1]),
                           'lsfr_bin': pd.Interval(args.lsfr_range[0], args.lsfr_range[1]),
                           'medlsm': median_lsm,
                           'medsfr': median_sfr,
                           'medlgm': median_lgm,
                           'medlgsm': np.nanmedian(lgsm),
                           'ledssfr9': np.nanmedian(np.log10(df_snap["sfr"]) - df_snap["log_star_mass"]) + 9,
                           'medgsfr9': np.nanmedian(np.log10(df_snap["sfr"]) - df_snap["log_gas_mass"]) + 9,
                           'medgssfr9': np.nanmedian(np.log10(df_snap["sfr"]) - lgsm) + 9,
                           'r0': pars[0],
                           'er0': epars[0],
                           'D0': pars[1],
                           'eD0': epars[1],
                           }, index=[0])
    
    outfile = f'{args.resfile_prefix}_inc_{args.inc_range[0]}_{args.inc_range[1]}.txt'
    df_out.to_csv(outfile, mode='a', sep='\t', header=not os.path.exists(outfile), index=None)

    return

# -----main code-----------------
if __name__ == '__main__':
    #	--------------------------	Read inputs	-------------------------------
    args = parse_args()

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
            execute_mode_indi(df_snap, args)
        elif (args.mode=='halo'):
            execute_mode_halo(df_snap, args)            
        elif (args.mode=='lsmzsfr'):
            execute_mode_lsmzsfr(df_snap, args)
        else:
            print("\nHmm...What mode is that again...?\n")

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))











































