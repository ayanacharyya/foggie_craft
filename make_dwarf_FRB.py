#!/usr/bin/env python3

"""
    Title :      make_dwarf_FRB
    Notes :      Compute DM contributions along various LoS from a SINGLE snapshot fo a dwarf galaxy
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run make_dwarf_FRB.py --system ayan_ssd --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 200 --res 0.5
                 run make_dwarf_FRB.py --system ayan_pleiades --foggie_dir /nobackupp19/aachary2/LowZRuns/ --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 200 --res 0.5
"""
from foggie_header import *
from make_3D_FRB_electron_density import plot_projection_diskrel, get_AM_vector, FITSImageData
from get_foggie_metallicity_profile import my_foggie_load, get_halo_coords, quant_dict
start_time = datetime.now()

# -----main code-----------------
if __name__ == '__main__':
    args = parse_args()
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    quant_arr = ['el_density', 'density']

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.fits_dir = args.output_dir + 'data/'
    Path(args.fits_dir).mkdir(parents=True, exist_ok=True)

    args.res = args.res_arr[0]
    args.res_text = f'_res{args.res:.1f}kpc'
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'

    fitsname = Path(args.fits_dir) / f'{args.output}_{args.halo}_FRB_{quant_dict[quant_arr[0]][0]}{args.upto_text}{args.res_text}.fits'

    # --------------getting stellar masses--------------
    input_filename = args.foggie_dir + 'LowZRunData.txt'
    df = pd.read_csv(input_filename, sep='\s+')

    # -----------making FRB------------------
    if not os.path.exists(fitsname) or args.clobber:
        # ----------reading in snapshot-----------
        halos_df_name = Path(args.code_path) / f'halo_infos/00{args.halo}/{args.run}/halo_c_v'
        halo_center = get_halo_coords(halos_df_name, args.output)
        snap_name = Path(args.foggie_dir) / f'halo_{int(args.halo):06d}' / args.run / args.output / args.output        
        
        ds, _ = my_foggie_load(snap_name, halo_center)
        args.current_redshift = ds.current_redshift
        args.current_time = ds.current_time.in_units('Gyr').tolist()

        sfr = -99 # dummy value for now
        log_star_mass = df['DD' + df['DD'].astype(str) == args.output]['log(M*/Msun)'].iloc[0]

        # ----------determining extent-----------------------
        if args.upto_kpc is not None:
            if args.docomoving:
                args.galrad = args.upto_kpc / (1 + args.current_redshift) / 0.695  # fit within a fixed comoving kpc h^-1, 0.695 is Hubble constant
            else:
                args.galrad = args.upto_kpc  # fit within a fixed physical kpc
        else:
            args.re = get_re_from_coldgas(args) if args.use_gasre else get_re_from_stars(ds, args)
            args.galrad = args.re * args.upto_re  # kpc

        # -------extract the required box------------
        box_center = ds.halo_center_kpc
        box_width = args.galrad * 2  # in kpc
        args.ncells = int(box_width / args.res)
        args.kpc_per_pix = 2 * args.galrad / args.ncells
        
        box_width_kpc = ds.arr(box_width, 'kpc')
        box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2.]

        # ----------------making projection plots-------------------
        norm_L = get_AM_vector(ds, radius_kpc=3., use_particles=False)

        # -------making and plotting the 3D FRBs--------------
        all_data = ds.arbitrary_grid(left_edge=box.left_edge, right_edge=box.right_edge, dims=[args.ncells, args.ncells, args.ncells])
        img_hdu_list = []

        for index, quant in enumerate(quant_arr):
            print(f'Making and plotting FRB for {quant} which is {index+1} out of {len(quant_arr)} quantities..')

            # --------making the 3D FRB------------
            FRB = all_data[('gas', quant_dict[quant][0])].in_units(quant_dict[quant][2]).astype(np.float32)

            # --------making the FITS ImageHDU---------------
            img_hdu = FITSImageData(FRB, ('gas', quant_dict[quant][1]))
            header = img_hdu[0].header
            for ind in range(3):
                header[f'CDELT{ind+1}'] = args.kpc_per_pix
                header[f'CUNIT{ind+1}'] = 'kpc'
                header[f'NORMAL_UNIT_VECTOR{ind+1}'] = norm_L[ind]

            header[f'SFR'] = sfr
            header[f'SFRUNIT'] = 'Msun/yr'

            header[f'LOG_MSTAR'] = log_star_mass
            header[f'MSTARUNIT'] = 'Msun'

            header[f'REDSHIFT'] = args.current_redshift

            img_hdu_list.append(img_hdu)

            # ------making the plots-----------
            fig_diskrel = plot_projection_diskrel(box, quant_dict[quant][0], box_width, norm_L, args, quant_label=quant_dict[quant][0], unit='Msun/pc**2' if quant == 'density' else quant_dict[quant][2], clim=[quant_dict[quant][3], quant_dict[quant][4]],  cmap=quant_dict[quant][6])

        # ------saving fits file------------------
        combined_img_hdu = FITSImageData.from_images(img_hdu_list)
        combined_img_hdu.writeto(fitsname, overwrite=args.clobber)
        print(f'Saved fits file as {fitsname}')
    else:
        print(f'Skipping snapshot {args.output} as {fitsname} already exists. Use --clobber to remake figure.')
        combined_img_hdu = fits.open(fitsname)

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
