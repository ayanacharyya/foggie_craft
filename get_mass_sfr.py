#!/usr/bin/env python3

"""

    Title :      get_mass_sfr
    Notes :      Derive stellar masses and SFRs for a list of FOGGIE snapshots
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    Jan 2026
    Examples :   run get_mass_sfr.py --system ayan_pleiades
"""
from header import *
start_time = datetime.now()

# -----------------------------------------------------------------------------
def get_sfr_df(args):
    '''
    Reads in the dataframe that contains SFR of each snapshot of a given halo
    Returns pandas dataframe
    '''
    sfr_filename = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/sfr'
    if os.path.exists(sfr_filename):
        print('Reading SFR history from', sfr_filename)
        sfr_df = pd.read_table(sfr_filename, names=('output', 'redshift', 'sfr'), comment='#', delim_whitespace=True)
    else:
        print('Did not find', sfr_filename, ', therefore will not include SFR')
        sfr_df = pd.DataFrame()

    return sfr_df

# -----main code-----------------
if __name__ == '__main__':
    args_tuple = parse_args('8508', 'RD0042')  # default simulation to work upon when comand line args not provided
    if type(args_tuple) is tuple: args, ds, refine_box = args_tuple # if the sim has already been loaded in, in order to compute the box center (via utils.pull_halo_center()), then no need to do it again
    else: args = args_tuple
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'

   # ---------initialising output dataframe-------------
    output_dfname = args.output_dir + 'data/mass_sfr_table.txt'
    df_out = pd.DataFrame(columns=['halo', 'snap', 'redshift', 'log_mass', 'sfr'])

    # ----------getting list of snapshots-----------
    halos = ['8508', '5036', '5016', '4123', '2392', '2878']
    redshift_list = [4, 3, 2, 1.5, 1, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]

    # ---------looping over halos-----------
    for thishalo in halos:
        print(f'Starting halo {thishalo}..')
        args.halo = thishalo
        sfr_df  = get_sfr_df(args) # reading SFR df
        
        # --------determining snapshots-----------
        df = pd.read_csv(args.code_dir + f'halo_infos/00{thishalo}/nref11c_nref9f/halo_cen_smoothed', sep=r'\s*\|\s*', engine='python')
        df = df.dropna(axis=1, how='all')[['snap', 'redshift']]
        output_list = []
        for redshift in redshift_list:
            idx = (df['redshift'] - redshift).abs().idxmin()
            output_list.append(df.loc[idx, 'snap'])

        # ---------looping over snapshots-----------
        for thisoutput in output_list:
            print(f'Starting snapshot {thishalo}:{thisoutput}..')
            args.output = thisoutput

            # -----------determining SFR amd stellar mass and redshift--------------------
            log_mstar = np.log10(get_disk_stellar_mass(args))
            try:
                sfr = sfr_df[sfr_df['output'] == args.output]['sfr'].values[0]
                redshift = sfr_df[sfr_df['output'] == args.output]['z'].values[0]
            except:
                sfr = -99
                redshift = -99

            # ------appending to dataframe----------
            df_out.loc[len(df_out)] = [thishalo, thisoutput, redshift, log_mstar, sfr]

    # -----------saving dataframe---------------
    df_out.to_csv(output_dfname, index=None, mode='a', header=not os.path.exists(output_dfname))
    print(f'Saved dataframe as {output_dfname}')
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds), args)
