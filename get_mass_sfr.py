#!/usr/bin/env python3

"""

    Title :      get_mass_sfr
    Notes :      Derive stellar masses and SFRs for a list of FOGGIE snapshots
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    Jan 2026
    Examples :   run get_mass_sfr.py --system ayan_pleiades --upto_kpc 200
"""
from foggie_header import *
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

# -----------------------------------------------------------------------------
def get_mass_profile(args):
    '''
    Function to get the mass profile for a stars and gas, for a given output
    '''
    mass_filename = args.code_path + 'halo_infos/00' + args.halo + '/' + args.run + '/masses_z-less-2.hdf5'

    if os.path.exists(mass_filename):
        print('Reading in', mass_filename)
        alldata = pd.read_hdf(mass_filename, key='all_data')
        thisdata = alldata[alldata['snapshot'] == args.output]

        if len(thisdata) == 0: # snapshot not found in masses less than z=2 file, so try greater than z=2 file
            mass_filename = mass_filename.replace('less', 'gtr')
            print('Could not find spanshot in previous file, now reading in', mass_filename)
            alldata = pd.read_hdf(mass_filename, key='all_data')
            thisdata = alldata[alldata['snapshot'] == args.output]

            if len(thisdata) == 0: # snapshot still not found in file
                print('Snapshot not found in either file. Returning bogus mass')
                return np.nan

        mass_profile = thisdata[thisdata['radius'] <= args.massrad]
        if len(mass_profile) == 0: # the smallest shell available in the mass profile is larger than the necessary radius within which we need the stellar mass
            print('Smallest shell avialable in mass profile is too small compared to args.massrad. Returning bogus mass')
            return np.nan
    else:
        print('File not found:', mass_filename)
        return np.nan

    return mass_profile

# -----------------------------------------------------------------------------
def get_masses_and_re(args, get_re_using='gas_HI_mass'):
    '''
    Function to determine the stellar, gas, and halo (total) mass, for a given snapshot, which is defined as the mass contained within args.massrad, which can either be a fixed absolute size in kpc OR = args.upto_re*Re
    Also determines the effective radius of stellar disk, based on the stellar mass profile
    Returns halo mass, gas mass, stellar mass, and stellar half-light radius
    '''

    mass_profile = get_mass_profile(args)
    mass_profile = mass_profile.sort_values('radius')

    mhalo = mass_profile['total_mass'].iloc[-1] # radius is in kpc, mass in Msun
    mgas = mass_profile['gas_mass'].iloc[-1] # radius is in kpc, mass in Msun
    mstar = mass_profile['stars_mass'].iloc[-1] # radius is in kpc, mass in Msun

    total_mass = mass_profile[get_re_using].iloc[-1]
    half_mass_radius = mass_profile[mass_profile[get_re_using] <= total_mass/2]['radius'].iloc[-1]

    return mhalo, mgas, mstar, half_mass_radius

# -----main code-----------------
if __name__ == '__main__':
    args_tuple = parse_args('8508', 'RD0042')  # default simulation to work upon when comand line args not provided
    if type(args_tuple) is tuple: args, ds, refine_box = args_tuple # if the sim has already been loaded in, in order to compute the box center (via utils.pull_halo_center()), then no need to do it again
    else: args = args_tuple
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'
    
   # ---------initialising output dataframe-------------
    output_dfname = args.output_dir + 'data/lsm_sfr_upto_disk.txt'
    df_out = pd.DataFrame(columns=['halo', 'snap', 'redshift', 're', 'disk_rad', 'log_star_mass', 'sfr', 'log_gas_mass', 'log_halo_mass'])

    # ----------getting list of snapshots-----------
    halos = ['8508', '5036', '5016', '4123', '2392', '2878']
    redshift_list = [5.0, 4.8, 4.6, 4.4, 4.2, 4.0,  3.8, 3.6, 3.4, 3.2, 3.0, 2.8, 2.6, 2.4, 2.2, 2.0, 1.8, 1.6, 1.4, 1.2, 1.0, 0.8, 0.6, 0.4, 0.2, 0.]

    # ---------looping over halos-----------
    for thishalo in halos:
        print(f'\nStarting halo {thishalo}..')
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
            print(f'\nStarting snapshot {thishalo}:{thisoutput}..')
            args.output = thisoutput

            # -----------determining SFR and redshift--------------------
            if args.output in sfr_df['output'].values:
                sfr = sfr_df[sfr_df['output'] == args.output]['sfr'].values[0]
                args.current_redshift = sfr_df[sfr_df['output'] == args.output]['redshift'].values[0]

                if args.docomoving: args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
                else: args.galrad = args.upto_kpc  # fit within a fixed physical kpc

                # ------determining extent for computing mass--------
                args.massrad = get_disk_rad(args)

                # ------determining stellar mass--------                
                mhalo, mgas, mstar, half_mass_radius = get_masses_and_re(args, get_re_using='gas_HI_mass')
                log_mhalo = np.log10(mhalo)
                log_mstar = np.log10(mstar)
                log_mgas = np.log10(mgas)

                # ------appending to dataframe----------
                df_out.loc[len(df_out)] = [thishalo, thisoutput, args.current_redshift, half_mass_radius, args.massrad, log_mstar, sfr, log_mgas, log_mhalo]
            else:
                print(f'Snapshot {args.output} is not in sfr_df, so filling dataframe with dummy values for this snapshot')
                df_out.loc[len(df_out)] = [thishalo, thisoutput, -99, -99, -99, -99, -99, -99, -99]

    # -----------saving dataframe---------------
    with open(output_dfname, 'w') as f:
        f.write('#') 
        df_out.to_csv(f, sep='\t', index=False)
    print(f'Saved dataframe as {output_dfname}')
    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
