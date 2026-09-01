#!/usr/bin/env python3

"""
    Title :      combine_foggie_vertical_metallicity_profile
    Notes :      Plot metallicity profile of a given list of filenames (FOGGIE snapshots)
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run combine_foggie_vertical_metallicity_profile.py --halo 2392,4123,5016 --upto_kpc 10
"""
from craft_header import *
from craft_utils import *
setup_plot_style()

start_time = datetime.now()

# ---------------------main code-----------------------------
if __name__ == '__main__':
    args = parse_args()  # default simulation to work upon when comand line args not provided
    if not args.keep: plt.close('all')
    args.fontfactor = 1.
    args.fontsize = 15
    if not 'pleiades' in args.system: args.output_dir = args.output_dir.replace('CRAFT', 'sharing')

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'metallicity_plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.data_dir = Path(args.output_dir) / 'metallicity_data'
    args.data_dir.mkdir(exist_ok=True, parents=True)
    
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'

    files = []
    for halo in args.halo_arr:
        thisfiles = glob.glob(str(args.data_dir) + f'/*_{halo}_binned_metaledge_on{args.upto_text}.txt')
        files.append(thisfiles)

    # ---------------setting up the plot----------------
    files = np.ravel(files)
    fig, ax = plt.subplots(figsize=(8, 6))


    # ---------------getting the metallicity profiles----------------
    z_minor_list = []
    z_major_list = []
    radius_grid = None
    minor_color = 'cornflowerblue'
    major_color = 'salmon'

    for index, file in enumerate(files):
        print(f'Processing file ({index + 1}/{len(files)}): {file}..')
        df = pd.read_csv(file, sep=r'\s+')
        
        # Grab and standardize radius grid from the first file
        if radius_grid is None:
            radius_grid = np.round(df['radius_kpc'].to_numpy(), decimals=5)
        
        z_minor_list.append(df['Z_minor'].to_numpy())
        z_major_list.append(df['Z_major'].to_numpy())

        # Plotting individual profiles for each halo
        ax.plot(radius_grid, df['Z_minor'], c=minor_color, alpha=0.5, lw=0.5)
        ax.plot(radius_grid, df['Z_major'], c=major_color, alpha=0.5, lw=0.5)

    z_minor_arr = np.array(z_minor_list)
    z_major_arr = np.array(z_major_list)

    df_summary = pd.DataFrame({
        'radius_kpc': radius_grid,
        
        'Z_minor_median': np.median(z_minor_arr, axis=0),
        'Z_minor_p16': np.percentile(z_minor_arr, 16, axis=0),
        'Z_minor_p84': np.percentile(z_minor_arr, 84, axis=0),
        'Z_minor_std': np.std(z_minor_arr, axis=0),
        
        'Z_major_median': np.median(z_major_arr, axis=0),
        'Z_major_p16': np.percentile(z_major_arr, 16, axis=0),
        'Z_major_p84': np.percentile(z_major_arr, 84, axis=0),
        'Z_major_std': np.std(z_major_arr, axis=0),
    })

    # Plotting the median profiles with shaded regions for percentiles
    ax.plot(radius_grid, df_summary['Z_minor_median'], c=minor_color, lw=2, label='Minor Axis')
    ax.fill_between(radius_grid, df_summary['Z_minor_p16'], df_summary['Z_minor_p84'], color=minor_color, alpha=0.3)
    ax.plot(radius_grid, df_summary['Z_major_median'], c=major_color, lw=2, label='Major Axis')
    ax.fill_between(radius_grid, df_summary['Z_major_p16'], df_summary['Z_major_p84'], color=major_color, alpha=0.3)

    ax = annotate_axes(ax, 'Radius (kpc)', r'Log metallicity (Z/Z$_\odot$)', args=args)
    ax.legend(fontsize=args.fontsize / args.fontfactor, loc='best', frameon=False)

    # ------------saving figure--------------
    figname = f'combined_binned_metaledge_on{args.upto_text}.png'
    save_fig(fig, Path(args.fig_dir), figname, args)
    
    # ------------saving txtfile--------------
    outfilename = args.data_dir / f'combined_binned_metaledge_on{args.upto_text}.csv'
    df_summary.to_csv(outfilename, index=False)
    print(f'Saved combined radial profile in {outfilename}')

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))

