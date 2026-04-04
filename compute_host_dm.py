#!/usr/bin/env python3
"""
    Title :      compute_host_dm
    Notes :      Estimate the observed DM for a host galaxy, given FRB/host properties and the existing FOGGIE models
    Output :     Plots as PDF
    Author :     Ayan Acharyya
    Started :    04-04-26
    Examples :   run compute_host_dm.py --use_sfr --input_cat frbcat0.txt
                 run compute_host_dm.py --input_cat frbcat0.txt
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
start_time = datetime.now()

# ------------------------------------------------------------------------------------------------
def make_latex_table(df_dmpars, outfilename, args, columns_to_publish=['id', 'lsm', 'impf', 'redshift', 'dm_16', 'dm_50', 'dm_84', 'dm_fit']):
    '''
    Convert the input dataframe into a latex table
    Saves latex table
    Returns latex dataframe
    '''    
    colnames_dict = {'id': 'FRB',
                     'lsm':r'\makecell{$\log(M_*/M_\odot$)\\range}', 
                     'impf':r'\makecell{Offset\\(kpc)}', 
                     'redshift':r'\makecell{Redshift}', 
                     'dm_16':r'\makecell{$DM_{16}$\\($pc\: cm^{-3}$)}',
                     'dm_50':r'\makecell{$DM_{50}$\\($pc\: cm^{-3}$)}',
                     'dm_84':r'\makecell{$DM_{84}$\\($pc\: cm^{-3}$)}',
                     'dm_fit':r'\makecell{$DM_{\rm fit}$\\($pc\: cm^{-3}$)}',
                     }

    columns_onedec = ['lsm', 'impf']
    columns_fourdec = ['redshift']

    df_latex = df_dmpars[columns_to_publish]
    df_mread = df_latex.copy()   
    
    for col in columns_onedec:
        df_latex[col] = df_latex[col].map('{:.1f}'.format)

    for col in columns_fourdec:
        df_latex[col] = df_latex[col].map('{:.4f}'.format)

    for col in (set(df_latex.columns) - set(np.hstack([columns_onedec, columns_fourdec, ['id']]))):
        df_latex[col] = df_latex[col].map('{:.0f}'.format)

    df_latex = df_latex.rename(columns=colnames_dict)

    df_mread.to_csv(str(outfilename).replace('.tex', '.txt'), index=None, sep='\t')
    df_latex.to_latex(outfilename, index=False, escape=False, column_format='l' * 1 + 'c' * (len(df_latex.columns) - 1))
    
    # -----------to insert lines between SFR groups------------
    insert_line_in_file('\\toprule\n', 1, outfilename) # to insert an additioal \toprule
    print(f'Saved latex table as {outfilename} and as .txt')

    return df_latex
 
# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')

    # -----------------read list of simulation outputs--------------
    filespecs =	args.data_dir / "lsm_sfr.txt"
    df_snap = pd.read_csv(filespecs, sep=r'\s+', engine='python')
    df_snap = df_snap.rename(columns={f'{df_snap.columns[0]}':f'{df_snap.columns[0][1:]}'})
    df_snap['log_sfr'] = np.log10(df_snap['sfr'])

    # ---------setting up filenames------------------
    if not os.path.exists(args.input_cat):
        args.input_cat = args.root_dir / args.input_cat

    sfr_text = '_with_sfr' if args.use_sfr else ''
    outfilename = args.input_cat.parent / f'{args.input_cat.stem}_with_dm{sfr_text}.txt'

    # ---------reading observation catalog------------------
    if not os.path.exists(outfilename) or args.clobber:
        print(f'Starting LoS DM calculations..')
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
        inverted_col_dict = {v: k for k, v in input_column_dict.items()}
        df_obs = pd.read_table(args.input_cat, sep=r'\s+')
        df_obs = df_obs.rename(columns = inverted_col_dict)
        print(f'Read in observed data from {args.input_cat}')
        #df_obs = df_obs[:5]

        # ------------declaring dummy columsn to be filled inside the loop---------------
        for col in ['sfr_low', 'sfr_med', 'sfr_up']:
            df_obs[col] = np.log10(df_obs[col])
        df_obs['dm_50'] = np.nan
        df_obs['dm_16'] = np.nan
        df_obs['dm_84'] = np.nan
        df_obs['dm_fit'] = np.nan

        # -------------looping over observed data------------------
        for index, obs in df_obs.iterrows():
            print(f'\nDoing object {obs["id"]} ({index + 1}/{len(df_obs)})...')
            # -----------------determine stellar mass range--------------------------
            lsm_range = [obs['lsm'] - obs_lsm_allowance, max(obs['lsm'] + obs_lsm_allowance, min_max_lsm)]

            # ----------------determine SFR range-------------------
            lsfr_range = [-10, 10]
            if args.use_sfr:
                if np.isfinite(obs['sfr_up']) and np.isfinite(obs['sfr_low']):
                    lsfr_range = [obs['sfr_low'], obs['sfr_up']]
                elif np.isfinite(obs['sfr_med']):
                    lsfr_range = [obs['sfr_med'] - obs_lsm_allowance, obs['sfr_med'] + obs_lsm_allowance]
                elif np.isfinite(obs['sfr_up']):
                    lsfr_range[1] = obs['sfr_up']

            # ------------find snapshots--------------------------
            df_sub = df_snap[(df_snap['log_star_mass'].between(lsm_range[0], lsm_range[1])) & 
                                (df_snap['log_sfr'].between(lsfr_range[0], lsfr_range[1]))].reset_index(drop=True)
            
            if len(df_sub) < 1:
                print(f'\n\nCAUTION: No simulated snapshots found for id {obs["id"]} in log SFR range {lsfr_range}, so skipping this object. Please re-run with increased obs_lsfr_allowance parameter in globalpars.py\n\n')
                continue

            print (f'\t\tFound {len(df_sub)} snapshots, within log mass range {lsm_range}, log sfr range {lsfr_range}')
            
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

            # --------------------reading in relevant LoS files of relevant snapshots-------------------
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

            # --------------estimate median DM contribiution--------------------
            stats = np.nanpercentile(combined_df['losdm'], (16, 50, 84))/2.0
            df_obs.at[index, 'dm_16'], df_obs.at[index, 'dm_50'], df_obs.at[index, 'dm_84'] = list(stats)


            # --------------------------perform fitting------------------------
            if args.use_sfr and np.isfinite(obs['sfr_med']):
                fit_D0 = 10 ** (lsm_lsfr_fit_D0[0] * (obs['lsm'] - 10) + lsm_lsfr_fit_D0[1] * obs['sfr_med'] + lsm_lsfr_fit_D0[2])     
                fit_r0 = 10 ** (lsm_lsfr_fit_r0[0] * (obs['lsm'] - 10) + lsm_lsfr_fit_r0[1] * obs['sfr_med'] + lsm_lsfr_fit_r0[2])       
            else:
                fit_D0 = 10.0 ** np.poly1d(lsm_fit_D0)(obs['lsm'] - 10)
                fit_r0 = 10.0 ** np.poly1d(lsm_fit_r0)(obs['lsm'] - 10)
            
            df_obs.at[index, 'dm_fit'] = 10 ** logradialexp3(obs['impf'], fit_r0, fit_D0)

        # ------------save output files---------------
        df_obs.to_csv(outfilename, sep='\t', index=None)
        print(f'Saved output as {outfilename}')
    else:
        print(f'Reading from existing {outfilename}; use --clobber to over-write')
    
    # -----------Read output file-----------------
    df_obs = pd.read_table(outfilename, sep='\t')
    texfilename =  args.fig_dir / Path(str(Path(outfilename).stem) + '.tex')
    df_latex = make_latex_table(df_obs, texfilename, args)

    # -------------DM vs DM plot------------------
    fig, ax = plt.subplots(figsize=(3, 2.8), layout='constrained')
    ax.errorbar(df_obs['dm_fit'], df_obs['dm_50'], yerr=[df_obs['dm_50'] - df_obs['dm_16'], df_obs['dm_84'] - df_obs['dm_50']], fmt='o', capsize=2, markersize=6)
    ax.plot([7, 120], [7, 120], ls='--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax = annotate_axes(ax, r'DM from scaling relation (pc cm$^{-3}$)', r'DM from LoS (pc cm$^{-3}$)', args=args, set_ticks=False)
    save_fig(fig, args.fig_dir, f'observed_DM{sfr_text}.pdf', args)

    # -------------DM distribution plot------------------
    fig, ax = plt.subplots(figsize=(3, 2.8), layout='constrained')
    delta_hist = 0.1
    ydata, bins, _ = ax.hist(np.log10(df_obs['dm_50']), bins=np.arange(1, 2 + delta_hist, delta_hist), color='b')

    bin_centers = (bins[1:] + bins[:-1]) / 2
    popt, pcov = curve_fit(gaussian, bin_centers, ydata)#, p0=initial_guess)
    #ax.plot(bin_centers, gaussian(bin_centers, *popt), color='r')
    ax.text(0.05, 0.95, r'$\mu$=' + f'{popt[1]:.2f}', color='k', transform=ax.transAxes, ha='left', va='top', fontsize=args.fontsize / args.fontfactor)
    ax.text(0.05, 0.85, r'$\sigma$=' + f'{popt[2]:.2f}', color='k', transform=ax.transAxes, ha='left', va='top', fontsize=args.fontsize / args.fontfactor)

    ax = annotate_axes(ax, r'DM from LoS (pc cm$^{-3}$)', 'Frequency', args=args, set_ticks=False)
    save_fig(fig, args.fig_dir, f'DM_distribution{sfr_text}.pdf', args)

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
