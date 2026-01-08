from __future__ import print_function

def get_run_loc_etc(args):
    if args.system == "ayan_local":
        foggie_dir = "/Users/acharyya/models/simulation_output/foggie/"
        output_path = "/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/FOGGIE_CRAFT/"
        code_path = "/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/"
    elif args.system == "ayan_hd":
        #foggie_dir = "/Volumes/Elements/foggieforayan/"
        foggie_dir = "/Volumes/Elements/acharyya_backup/models/simulation_output/foggie/"
        output_path = "/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/FOGGIE_CRAFT/"
        code_path = "/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/"
    elif args.system == "ayan_pleiades":
        foggie_dir = "/nobackup/mpeeples/" if args.foggie_dir is None else args.foggie_dir
        output_path = "/nobackupp19/aachary2/foggie_craft/"
        code_path = "/nobackupp19/aachary2/foggie_craft/foggie/foggie/"

    if not args.pwd:
        if args.run == "natural":
            runname = "nref11n"
        elif args.run == "nref10f" or args.run == "nref11n_nref10f":
            runname = "nref11n_nref10f"
        elif args.run == "nref11c_nref9f" or args.run == "nref11c":
            runname = "nref11c_nref9f"
        elif ('feedback' in args.run) and (not 'track' in args.run) and (not args.forcepath):
            runname = "nref11c_nref9f"
        else:
            runname = args.run

        trackname = code_path + "halo_tracks/00" + args.halo + "/nref11n_selfshield_15/halo_track_200kpc_nref9"
        infofile = code_path  + "halo_infos/00" + args.halo + "/" + runname + "/halo_info"
        haloname = "halo_00" + args.halo + "_" + runname

        if ('feedback' in args.run) and (not 'track' in args.run) and (not args.forcepath):
            run_loc = args.run + "/"
        else:
            run_loc = "halo_00" + args.halo + "/" + runname + "/"

        output_dir = output_path #+ "plots_" + run_loc
        spectra_dir = output_dir + "spectra/"

        if args.system=='cassiopeia' or args.system=='pleiades_cassi':
            output_dir = output_path
            runname = args.run
            run_loc = "halo_00" + args.halo + "/" + runname + "/"

    if args.pwd:
        print('using pwd args')
        foggie_dir = '.'
        output_path = '.'
        output_dir = './'
        run_loc = '.'
        code_path = '.'
        trackname = 'halo_track'
        haloname = 'halo'
        infofile = 'halo_info'
        spectra_dir = '.'

    return foggie_dir, output_dir, run_loc, code_path, trackname, haloname, spectra_dir, infofile
