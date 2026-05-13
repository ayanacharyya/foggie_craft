#!/usr/bin/env python3

"""
    Title :      get_foggie_vertical_metallicity_profile
    Notes :      Plot metallicity profile of a given list of filenames (FOGGIE snapshots)
    Output :     Pandas dataframe
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run get_foggie_vertical_metallicity_profile.py --system ayan_ssd --halo 8508 --run nref11c_nref9f --output RD0027 --upto_kpc 20
                 run get_foggie_vertical_metallicity_profile.py --system ayan_ssd --halo 2392 --run 2162_39 --output DD1360 --upto_kpc 20
"""
from craft_header import *
from craft_utils import *
setup_plot_style()
from foggie_header import *
from make_3D_FRB_electron_density import get_AM_vector
from get_foggie_metallicity_profile import get_halo_coords, my_foggie_load, get_projection_frb, quant_dict
start_time = datetime.now()

# ---------------------main code-----------------------------
if __name__ == '__main__':
    args = parse_args()  # default simulation to work upon when comand line args not provided
    if not args.keep: plt.close('all')
    args.fontfactor = 1.2
    if not 'pleiades' in args.system: args.output_dir = args.output_dir.replace('CRAFT', 'sharing')

    # -----------determining directories----------------
    args.fig_dir = args.output_dir + 'metallicity_plots/'
    Path(args.fig_dir).mkdir(parents=True, exist_ok=True)
    
    args.data_dir = Path(args.output_dir) / 'metallicity_data'
    args.data_dir.mkdir(exist_ok=True, parents=True)
    
    args.upto_text = '_upto%.1Fckpchinv' % args.upto_kpc if args.docomoving and args.upto_kpc is not None else '_upto%.1Fkpc' % args.upto_kpc if args.upto_kpc is not None else f'_upto{args.re:.1f}re'
    fitsname = args.data_dir / f'{args.halo}_{args.run}_{args.output}_projected_met_map{args.upto_text}.fits'

    # ---------------getting the metallicity profile----------------
    if not os.path.exists(fitsname) or args.clobber:
        print(f']\n{fitsname} does not exist. Creating afresh..')

        # ----------reading in snapshot-----------
        halos_df_name = Path(args.code_path) / f'halo_infos/00{args.halo}/{args.run}/halo_c_v'
        halo_center = get_halo_coords(halos_df_name, args.output)
        snap_name = Path(args.foggie_dir) / f'halo_{int(args.halo):06d}' / args.run / args.output / args.output        
        
        ds, _ = my_foggie_load(snap_name, halo_center)
        args.current_redshift = ds.current_redshift
        args.current_time = ds.current_time.in_units('Gyr').tolist()

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
        box_width_kpc = ds.arr(box_width, 'kpc')
        box = ds.r[box_center[0] - box_width_kpc / 2.: box_center[0] + box_width_kpc / 2., box_center[1] - box_width_kpc / 2.: box_center[1] + box_width_kpc / 2., box_center[2] - box_width_kpc / 2.: box_center[2] + box_width_kpc / 2.]

        # ----------------setting up FITS file-------------------
        norm_L = get_AM_vector(ds, radius_kpc=3., use_particles=False)
        quant_arr = ['metal', 'density']
        img_hdu_list = []
        primary_hdu = fits.PrimaryHDU()
        img_hdu_list.append(primary_hdu)

        plot_width = box_width/1.44
        ncell_buff = 800 # this buffer size is purely for visualising and saving the projected map
        kpc_per_pix_proj = plot_width / ncell_buff

        for quant in quant_arr:
            # ----------------making projection plots-------------------
            fig_diskrel, data_faceon, data_edgeon = get_projection_frb(box, quant_dict[quant][0], plot_width, norm_L, args, 
                                                    quant_label=quant_dict[quant][0], 
                                                    unit='Msun/pc**2' if quant == 'density' else quant_dict[quant][2], 
                                                    clim=[quant_dict[quant][3], quant_dict[quant][4]] if quant_dict[quant][3] is not None else None,  
                                                    cmap=quant_dict[quant][6], 
                                                    takelog=quant_dict[quant][7], 
                                                    clabel=quant_dict[quant][9],
                                                    return_frb=True,
                                                    do_plot=True,
                                                    ncell_buff = ncell_buff,
                                                    )
        
            # --------making the FITS ImageHDU for 2D projected FRB: face on---------------
            hdu_faceon = fits.ImageHDU(data=data_faceon.astype(np.float32))
            header = hdu_faceon.header
            header['EXTNAME'] = f'{quant_dict[quant][1]} FACE ON PROJ'
            header['FIELD'] = quant_dict[quant][1]
            header['BUNIT'] = quant_dict[quant][8]
            header['WIDTH'] = (plot_width, 'kpc')
            for ind in range(2):
                header[f'CDELT{ind+1}'] = kpc_per_pix_proj
                header[f'CUNIT{ind+1}'] = 'kpc'
            img_hdu_list.append(hdu_faceon)

            # --------making the FITS ImageHDU for 2D projected FRB: face on---------------
            hdu_edgeon = fits.ImageHDU(data=data_edgeon.astype(np.float32))
            header = hdu_edgeon.header
            header['EXTNAME'] = f'{quant_dict[quant][1]} EDGE ON PROJ'
            header['FIELD'] = quant_dict[quant][1]
            header['BUNIT'] = quant_dict[quant][8]
            header['WIDTH'] = (plot_width, 'kpc')
            for ind in range(2):
                header[f'CDELT{ind+1}'] = kpc_per_pix_proj
                header[f'CUNIT{ind+1}'] = 'kpc'
            img_hdu_list.append(hdu_edgeon)

        # -----------------saving the projections as fits file------------
        hdul = fits.HDUList(img_hdu_list)
        hdul.writeto(fitsname, overwrite=True)
        print(f'Successfully saved {len(hdul)-1} extensions to {fitsname}"')

    else:
        print(f'\nReading from existing file {fitsname}')
        hdul = fits.open(fitsname)

    # --------preparing cone bins------------
    quant = 'metal'
    data = hdul[f'{quant_dict[quant][1].upper()} EDGE ON PROJ'].data
    
    # ----------------making profile plots-------------------


    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
