#!/usr/bin/env python3
"""
    Title :      compute_host_dm
    Notes :      Estimate the observed DM for a host galaxy, given FRB/host properties and the existing FOGGIE models
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    04-04-26
    Examples :   run compute_host_dm.py --use_sfr --input_cat frbcat0.txt --plot_compare --plot_dist --plot_scaling
                 run compute_host_dm.py --input_cat frbcat0.txt --plot_all
                 run compute_host_dm.py --input_cat fgcat0.txt
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
start_time = datetime.now()

# --------------------------------------------------------------------------------------------------------------------
def get_SFMS_Popesso23(log_mass, redshift):
    '''
    Computes an empirical SFMS based on Popesso+23 (https://arxiv.org/abs/2203.10487) Eq 10, for given redshift
    Then returns a two-part log SFR array based on an input minimum and maximum log mass; the two part are based on the lower limit of this empirical relation and any extrapolation below it
    Returns two tuples: (log_mass1, log_SFR1) and (log_mass2, log_SFR2)
    '''
    a0, a1, b0, b1, b2 = ufloat(0.2,0.02), ufloat(-0.034, 0.002), ufloat(-26.134,0.015), ufloat(4.722, 0.012), ufloat(-0.1925, 0.0011)  # Table 2, Eq 10
    log_mass_low_lim = 8.7 # lower limit of mass they fitted up to

    age_at_z = cosmo.age(redshift).value # Gyr

    if log_mass < log_mass_low_lim: 
        log_SFR = (a1 * age_at_z + b1) * log_mass + b2 * (log_mass) ** 2 + b0 + a0 * age_at_z
    else: 
        log_SFR = (a1 * age_at_z + b1) * log_mass + b2 * (log_mass) ** 2 + b0 + a0 * age_at_z

    return log_SFR.n

# ------------------------------------------------------------------------------------------------
def read_foggie_catalog(filename):
    '''
    Function to read in the input FOGGIE catalog,
    Rename its columns if needed
    Returns pandas dataframe
    '''
    if filename.suffix == '.txt':
        df = pd.read_csv(filename, sep=r'\s+', engine='python')
    elif filename.suffix == '.csv':
        df = pd.read_csv(filename)

    df = df.drop_duplicates(subset=['halo', 'snap'], keep='last')

    if df.columns[0][0] == '#': # if the first column name ha a '#' then remove it
        df = df.rename(columns={f'{df.columns[0]}':f'{df.columns[0][1:]}'})

    if 'sfr' in df and 'log_sfr' not in df:
        df['log_sfr'] = np.log10(df['sfr'])

    df = df.rename(columns={'log_star_mass_from_snap': 'log_star_mass', 'log_gas_mass_from_profile': 'log_gas_mass'})
    df = df[(df['redshift'].between(args.z_range[0], args.z_range[1], inclusive='left'))].reset_index(drop=True)
    df['log_ssfr'] = df['log_sfr'] - df['log_star_mass']

    print (f'\t\tFound {len(df)} snapshots, within {args.z_range}')

    return df

# ------------------------------------------------------------------------------------------------
def read_obs_catalog(filename, args, input_column_dict=None, add_columns=['dm_16', 'dm_50', 'dm_84', 'dm_fit1d', 'dm_fit2d']):
    '''
    Function to read in the input observed catalog,
    Rename its columns as per input_column_dict, and
    Add new columns (add_columns) and set them to nan, to be filled in later
    Returns pandas dataframe
    '''
    if not os.path.exists(filename):
        filename = args.root_dir / filename
    
    if os.path.exists(filename):
        print(f'Reading in observed data from {filename}')
        df = pd.read_table(filename, sep=r'\s+', comment='#')

        if input_column_dict is not None:
            inverted_col_dict = {v: k for k, v in input_column_dict.items()}
            df = df.rename(columns = inverted_col_dict)

        # ------------declaring dummy columsn to be filled inside the loop---------------
        df = df.rename(columns={'SFR_low':'sfr_low', 'SFR_med':'sfr_med', 'SFR_up':'sfr_up', 'log(M_*)':'log_mass'})
        for col in ['sfr_low', 'sfr_med', 'sfr_up']:
            df[col] = np.log10(df[col])
        for col in add_columns:
            df[col] = np.nan
    else:
        print(f'Could not find {filename}..')
        df = pd.DataFrame()

    return df, filename

# ------------------------------------------------------------------------------------------------
def get_param_ranges(obs):
    '''
    Function to determine the ranges of SFR, stellar mass, inclination and impact factor in which to select snapshots and LoS
    Returns lsm_range, lsfr_range, inc_range, impf_range, each being a tuple of length 2
    '''
    # -----------------determine stellar mass range--------------------------
    lsm_range = [obs['lsm'] - obs_lsm_allowance, max(obs['lsm'] + obs_lsm_allowance, min_max_lsm)]

    # ----------------determine SFR range-------------------
    # lsfr_range = [-10, 10]
    # if np.isfinite(obs['sfr_up']) and np.isfinite(obs['sfr_low']):
    #     lsfr_range = [obs['sfr_low'], obs['sfr_up']]
    # elif np.isfinite(obs['sfr_med']):
    #     lsfr_range = [obs['sfr_med'] - obs_lsm_allowance, obs['sfr_med'] + obs_lsm_allowance]
    # elif np.isfinite(obs['sfr_up']):
    #     lsfr_range[1] = obs['sfr_up']
    lsfr_range = [obs['sfr_to_use'] - obs_lsm_allowance, obs['sfr_to_use'] + obs_lsm_allowance]
    
    # -----------------------determine inclination range----------------------
    if 'inc' in obs and np.isfinite(obs['inc']):
        inc_range = [obs['inc'] - obs_inc_allowance, obs['inc'] + obs_inc_allowance]
    else:
        inc_range = [0, 90]
    
    # -----------------------determine impact parameter range----------------------
    if (obs['impf_low'] + obs['impf_up']) >= obs['impf'] * obs_impf_frac_allow:
        impf_range = [obs['impf'] - obs['impf_low'], obs['impf'] + obs['impf_up']]
    else:
        impf_range = [obs['impf'] * (1 - obs_impf_frac_allow/2), obs['impf'] * (1 + obs_impf_frac_allow/2)]

    return lsm_range, lsfr_range, inc_range, impf_range

# ------------------------------------------------------------------------------------------------
def find_los_in_range(df_snap, inc_range, impf_range, args):
    '''
    Function to find and combine LoS within given inc_range and impf_range
    Returns combined dataframe containing all such LoS
    '''
    combined_df = pd.DataFrame()
    
    for i, snap in df_snap.iterrows():
        thisfile = args.los_dir / f'{snap["snap"]}_{snap["halo"]}_FRB_El_number_density_upto{args.rangekpc}kpc_res{args.reskpc}kpc_150.npy'
        if os.path.exists(thisfile):
            dm_arr	= np.load(thisfile)
        else:
            continue
        this_df = pd.DataFrame(dm_arr, columns=['inc', 'impf', 'distmaj', 'losdm'])
        this_df = this_df[(this_df['inc'].between(inc_range[0], inc_range[1])) & (this_df['impf'].between(impf_range[0], impf_range[1]))]
        combined_df = pd.concat([combined_df, this_df], ignore_index=True)
        
    print(f"Total number of LoS = {len(combined_df)}, within inclination range {inc_range}, impact factor range {impf_range}")	

    return combined_df

# ------------------------------------------------------------------------------------------------
def get_fit_r0_d0(obs, args):
    '''
    Function determine r0 and D0 from scaling relation fit, based on fit coeffciients that are read in from the *multifit_params.txt file
    Returns fit_D0, fit_r0
    '''
    paramfilename = f'{args.fig_dir}/{Path(args.resfile_prefix).stem}_z_{args.z_range[0]}_{args.z_range[1]}_DM0_r0_vs_lsm_inc_{args.inc_range[0]}_{args.inc_range[1]}_multifit_params.txt'
    results = np.loadtxt(paramfilename)
    lsm_lsfr_fit_D0 = results[0] # D0 vs SFR (100Myr) and mass
    lsm_lsfr_fit_r0 = results[2] # r0 vs SFR (100Myr) and mass

    lsm_fit_D0 = results[4][:-1] # D0 vs SFR (100Myr)
    lsm_fit_r0 = results[6][:-1] # r0 vs SFR (100Myr)

    fit2d_D0 = 10 ** (lsm_lsfr_fit_D0[0] * obs['sfr_to_use'] + lsm_lsfr_fit_D0[1] * (obs['lsm'] - 10) + lsm_lsfr_fit_D0[2])     
    fit2d_r0 = 10 ** (lsm_lsfr_fit_r0[0] * obs['sfr_to_use'] + lsm_lsfr_fit_r0[1] * (obs['lsm'] - 10) + lsm_lsfr_fit_r0[2])       

    fit1d_D0 = 10 ** np.poly1d(lsm_fit_D0)(obs['sfr_to_use'])
    fit1d_r0 = 10 ** np.poly1d(lsm_fit_r0)(obs['sfr_to_use'])

    return fit1d_D0, fit1d_r0, fit2d_D0, fit2d_r0

# ------------------------------------------------------------------------------------------------
def make_latex_table(df, outfilename, args, columns_to_publish=['id', 'lsm', 'sfr_med', 'impf', 'impf_low', 'impf_up', 'redshift', 'dm_16', 'dm_50', 'dm_84', 'dm_fit1d', 'dm_fit2d']):
    '''
    Convert the input dataframe into a latex table
    Saves latex table
    Returns latex dataframe
    '''    
    colnames_dict = {'id': 'FRB',
                     'lsm':r'\makecell{$\log(M_*/M_\odot$)}', 
                     'sfr_med':r'\makecell{$\log$(SFR /M$_\odot$ yr$^{\rm -1}$)}', 
                     'impf':r'\makecell{Offset\\(kpc)}', 
                     'redshift':r'\makecell{Redshift}', 
                     'dm':r'\makecell{DM$_{\rm los}$\\($pc\: cm^{-3}$)}',
                     'dm_fit1d':r'\makecell{DM$_{\rm fit, 1D}$\\($pc\: cm^{-3}$)}',
                     'dm_fit2d':r'\makecell{DM$_{\rm fit, 2D}$\\($pc\: cm^{-3}$)}',
                     }

    columns_with_err = ['impf', 'dm']
    columns_onedec = ['lsm']
    columns_threedec = ['redshift']
    columns_with_nans = ['sfr_med', 'dm_fit1d', 'dm_fit2d']

    df_latex = df[columns_to_publish]
    df_latex = df_latex.rename(columns={'dm_50':'dm', 'dm_16': 'dm_low', 'dm_84':'dm_up'})
    df_mread = df_latex.copy()   
    
    for col in columns_with_err:
        df_latex[col] = df_latex.apply(lambda row: rf'{row[col] :.1f} $^{{+{(row[col+"_up"] - row[col]) :.1f}}}_{{-{(row[col] - row[col+"_low"]) :.1f}}}$' if row[col+"_up"] > 0.1 and not np.isnan(row[col]) else rf'{row[col] :.1f}' if not np.isnan(row[col]) else '-', axis=1)
        df_latex.drop(columns=[col + '_low', col + '_up'], inplace=True)

    for col in columns_with_nans:
        df_latex[col] = df_latex.apply(lambda row: rf'{row[col] :.1f}' if not np.isnan(row[col]) else '-', axis=1)

    for col in columns_onedec:
        df_latex[col] = df_latex[col].map('{:.1f}'.format)

    for col in columns_threedec:
        df_latex[col] = df_latex[col].map('{:.3f}'.format)

    for col in (set(df_latex.columns) - set(np.hstack([columns_with_err, columns_onedec, columns_threedec, columns_with_nans, ['id']]))):
        df_latex[col] = df_latex[col].map('{:.0f}'.format)

    df_latex = df_latex.rename(columns=colnames_dict)

    df_mread.to_csv(str(outfilename).replace('.tex', '.txt'), index=None, sep='\t')
    df_latex.to_latex(outfilename, index=False, escape=False, column_format='l' * 1 + 'c' * (len(df_latex.columns) - 1))
    
    # -----------to insert lines between SFR groups------------
    insert_line_in_file('\\toprule\n', 1, outfilename) # to insert an additioal \toprule
    print(f'Saved latex table as {outfilename} and as .txt')

    return df_latex
 
# ------------------------------------------------------------------------------------------------
def plot_dm_comparison(df, args):
    '''
    Function to plot DM measured from LoS vs DM measured from scaling relations, from a given dataframe
    Saves the plot
    Returns figure handle
    '''
    fig, ax = plt.subplots(figsize=(3, 2.8), layout='constrained')
    ax.errorbar(df['dm_fit1d'], df['dm_50'], yerr=[df['dm_50'] - df['dm_16'], df['dm_84'] - df['dm_50']], fmt='o', capsize=2, markersize=6)
    ax.errorbar(df['dm_fit2d'], df['dm_50'], yerr=[df['dm_50'] - df['dm_16'], df['dm_84'] - df['dm_50']], fmt='s', capsize=2, markersize=6)
    ax.plot([7, 120], [7, 120], ls='--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax = annotate_axes(ax, r'DM from scaling relation (pc cm$^{-3}$)', r'DM from LoS (pc cm$^{-3}$)', args=args, set_ticks=False)
    save_fig(fig, args.fig_dir, f'DM_comparison.pdf', args)

    return fig

# ------------------------------------------------------------------------------------------------
def plot_dm_distribution(df, args):
    '''
    Function to plot distribution of DM measured from LoS, from a given dataframe
    Saves the plot
    Returns figure handle
    '''
    fig, ax = plt.subplots(figsize=(3, 2.8), layout='constrained')
    delta_hist = 0.1
    ydata, bins, _ = ax.hist(np.log10(df['dm_fit1d']), bins=np.arange(1, 2 + delta_hist, delta_hist), color='goldenrod', alpha=0.5)
    ydata, bins, _ = ax.hist(np.log10(df['dm_fit2d']), bins=np.arange(1, 2 + delta_hist, delta_hist), color='salmon', alpha=0.5)
    ydata, bins, _ = ax.hist(np.log10(df['dm_50']), bins=np.arange(1, 2 + delta_hist, delta_hist), color='cornflowerblue', alpha=0.5)

    bin_centers = (bins[1:] + bins[:-1]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, ydata)
    #ax.plot(bin_centers, gaussian(bin_centers, *popt), color='r') # to plot the fited distribution
    ax.text(0.05, 0.95, r'$\mu$=' + f'{popt[1]:.2f}', color='k', transform=ax.transAxes, ha='left', va='top', fontsize=args.fontsize / args.fontfactor)
    ax.text(0.05, 0.85, r'$\sigma$=' + f'{popt[2]:.2f}', color='k', transform=ax.transAxes, ha='left', va='top', fontsize=args.fontsize / args.fontfactor)

    ax = annotate_axes(ax, r'$\log$ (DM /pc cm$^{-3}$)', 'Number', args=args, set_ticks=False)
    save_fig(fig, args.fig_dir, f'DM_distribution.pdf', args)
    
    return fig

# ------------------------------------------------------------------------------------------------
def plot_dm_scaling(df, args):
    '''
    Function to plot scaling relations of DM measured from LoS vs various host properties, from a given dataframe
    Saves the plot
    Returns figure handle
    '''
    xcols = ['lsm', 'sfr_to_use', 'impf']
    ycol = 'dm_fit1d'
    label_dict = {'lsm':r'$\log(M_*/M_\odot$)', 'sfr_to_use':r'$\log$ SFR ($M_\odot$/yr)', 'impf':'Impact factor (kpc)'}

    df_high = df[df['lsm'] > 10.5]
    df_low = df[df['lsm'] <= 10.5]

    # --------------looping over each panel-----------
    fig, axes = plt.subplots(1, len(xcols), figsize=(9, 2.8), layout='constrained')
    for index, xcol in enumerate(xcols):
        #axes[index].errorbar(df[xcol], df['dm_50'], yerr=[df['dm_50'] - df['dm_16'], df['dm_84'] - df['dm_50']], fmt='o', capsize=2, markersize=6, markeredgecolor='k', markerfacecolor='cornflowerblue')
        axes[index].plot(df_high[xcol], df_high[ycol], 'bo', label=r'$\log(M_*/M_\odot$) > 10.5')
        axes[index].plot(df_low[xcol], df_low[ycol], 'bo', fillstyle='none', label=r'$\log(M_*/M_\odot$) < 10.5')

        axes[index].set_yscale('log')
        axes[index].set_yticks([1, 3, 10, 30], [1, 3, 10, 30])
        #if 'log' not in label_dict[xcol]: axes[index].set_xscale('log')
        axes[index] = annotate_axes(axes[index], label_dict[xcol], r'DM$_{\rm host}$ (pc cm$^{-3}$)', args=args, set_ticks=False, hide_yaxis=index)
        if index == 1:
            axes[index].legend(loc='upper left', fontsize=args.fontsize)
    
    save_fig(fig, args.fig_dir, f'DM_scaling.pdf', args)
    
    return fig

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    input_column_dict = {'id': 'FRB', 
                        'lsm': 'log(M_*)', 
                        'sfr_low': 'SFR_low', 
                        'sfr_med': 'SFR_med', 
                        'sfr_up': 'SFR_up', 
                        'impf': 'Offset', 
                        'impf_low': 'lower', 
                        'impf_up': 'higher', 
                        'redshift': 'z',
                        'inc': 'inc',
    }                                           # the user will need to modify this dict

    # ---------reading input catalogs------------------
    df_snap = read_foggie_catalog(args.data_dir / "lsm_sfr_masses_upto_disk.csv")
    add_columns = ['sfr_to_use', 'dm_16', 'dm_50', 'dm_84', 'dm_fit1d', 'dm_fit2d']
    df_obs, args.input_cat = read_obs_catalog(args.input_cat, args, input_column_dict=input_column_dict, add_columns=add_columns)

    # ---------setting output filenames------------------
    outfilename = args.input_cat.parent / f'{args.input_cat.stem}_with_dm.csv'

    # ---------making output catalog with DMs------------------
    if not os.path.exists(outfilename) or args.clobber:
        print(f'Starting LoS DM calculations..')
        
        # -------------looping over observed data------------------
        for index, obs in df_obs.iterrows():
            print(f'\nDoing object {obs["id"]} ({index + 1}/{len(df_obs)})...')
            # ---------derive SFR from SFMS if it does not exist----------------
            if np.isnan(obs['sfr_med']):
                obs['sfr_to_use'] = get_SFMS_Popesso23(obs['lsm'], obs['redshift'])
            else:
                obs['sfr_to_use'] = obs['sfr_med']
            df_obs.at[index, 'sfr_to_use'] = obs['sfr_to_use']

            # --------------------------derive fitted DM------------------------
            fit1d_D0, fit1d_r0, fit2d_D0, fit2d_r0 = get_fit_r0_d0(obs, args)
            df_obs.at[index, 'dm_fit1d'] = schechter(obs['impf'], fit1d_r0, fit1d_D0) # update dataframe
            df_obs.at[index, 'dm_fit2d'] = schechter(obs['impf'], fit2d_r0, fit2d_D0) # update dataframe

            # ------------find snapshots--------------------------
            lsm_range, lsfr_range, inc_range, impf_range = get_param_ranges(obs)
            df_sub = df_snap[(df_snap['log_star_mass'].between(lsm_range[0], lsm_range[1])) & (df_snap['log_sfr'].between(lsfr_range[0], lsfr_range[1]))].reset_index(drop=True)
            
            if len(df_sub) < 1:
                print(f'\n\nCAUTION: No simulated snapshots found for id {obs["id"]} in log SFR range {lsfr_range}, so skipping this object. Please re-run with increased obs_lsfr_allowance parameter in globalpars.py\n\n')
                continue
            print (f'\t\tFound {len(df_sub)} snapshots, within log mass range {lsm_range}, log sfr range {lsfr_range}')
            
            # --------------------reading in relevant LoS files of relevant snapshots-------------------
            combined_df = find_los_in_range(df_snap, inc_range, impf_range, args)
            stats = np.nanpercentile(combined_df['losdm'], (16, 50, 84))/2.0
            df_obs.at[index, 'dm_16'], df_obs.at[index, 'dm_50'], df_obs.at[index, 'dm_84'] = list(stats) # update dataframe


        # ------------save output files---------------
        df_obs.to_csv(outfilename, index=None)
        print(f'Saved output as {outfilename}')
    else:
        print(f'Reading output file from existing {outfilename}; use --clobber to over-write')
    
    # -----------Read output file-----------------
    df_output = pd.read_csv(outfilename)
    texfilename =  args.fig_dir / Path(str(Path(outfilename).stem) + '.tex')
    df_latex = make_latex_table(df_output, texfilename, args)

    # -------------making plots------------------
    if args.plot_compare or args.plot_all:
        fig_comp = plot_dm_comparison(df_output, args)
    if args.plot_distribution or args.plot_all:
        fig_dist = plot_dm_distribution(df_output, args)
    if args.plot_scaling or args.plot_all:
        fig_scaling = plot_dm_scaling(df_output, args)

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
