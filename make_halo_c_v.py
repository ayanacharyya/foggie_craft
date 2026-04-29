#!/usr/bin/env python3

"""
    Title :      make_halo_c_v
    Notes :      Reads a given file with info about dwarf runs and converts them into the format of halo_c_v files
    Author :     Ayan Acharyya
    Started :    29-04-2026
    Examples :   run make_halo_c_v.py --system ayan_ssd
                 run make_halo_c_v.py --system ayan_pleiades --foggie_dir /nobackupp19/aachary2/LowZRuns/
"""
from foggie_header import *
start_time = datetime.now()

# -----main code-----------------
if __name__ == '__main__':
    args_tuple = parse_args('8508', 'RD0042')  # default simulation to work upon when comand line args not provided
    if type(args_tuple) is tuple: args, ds, refine_box = args_tuple # if the sim has already been loaded in, in order to compute the box center (via utils.pull_halo_center()), then no need to do it again
    else: args = args_tuple
    if not args.keep: plt.close('all')
    if args.system == "ayan_pleiades": args.code_dir = '/nobackupp19/aachary2/ayan_codes/foggie/foggie/'
    else: args.code_dir = '/Users/acharyya/Work/astro/ayan_codes/foggie/foggie/'

    input_filename = args.foggie_dir + 'LowZRunData.txt'
    base_dir = args.code_dir + '/halo_infos/'

    halo_dict = {   '2392'  :  'Hurricane' ,
                    '2878'  :  'Cyclone' ,
                    '4123'  :  'Blizzard' ,
                    '5016'  :  'Squall' ,
                    '5036'  :  'Maelstrom' ,
                    '8508'  :  'Tempest',
                    }

    # Create a reverse map (Name -> 4-digit ID) to use in folder paths
    sim_to_id = {v: k for k, v in halo_dict.items() if len(k) == 4}

    # 2. Read the LowZRunData.txt file 
    # Using sep='\s+' to handle the irregular spacing in the .txt file
    df = pd.read_csv(input_filename, sep='\s+')

    # 3. Iterate through each row of the simulation data
    for index, row in df.iterrows():
        sim_name = row['Sim']
        
        # Determine the halo_id from the dictionary
        if sim_name in sim_to_id:
            halo_id = sim_to_id[sim_name]
        else:
            print(f"Warning: {sim_name} not found in halo_dict. Skipping...")
            continue
        
        snap_id = f"DD{int(row['DD']):04d}"
        run_id = f"{row['SatID']}"
        target_path = os.path.join(base_dir, f"00{halo_id}", run_id)
        
        if not os.path.exists(target_path):
            os.makedirs(target_path)
            print(f"Created directory: {target_path}")

        output_file = os.path.join(target_path, 'halo_c_v')
        
        # Header matching the template halo_c_v file
        header = "|            redshift |   name |               time |                x_c |                y_c |                z_c |                   v_x |                 v_y |                 v_z |"
        
        # Data values: Use Cen_x/y/z from table, others are -99 as requested
        # We use string formatting to match the fixed-width style of the original file 
        data_row = (f"| {-99.0:18.7f} | {snap_id:>6} | {-99.0:18.7f} | "
                    f"{row['Cen_x']:18.7f} | {row['Cen_y']:18.7f} | {row['Cen_z']:18.7f} | "
                    f"{-99.0:21.7f} | {-99.0:19.7f} | {-99.0:19.7f} |")
        
        if not os.path.exists(output_file):
            with open(output_file, 'w') as f:
                f.write(header + "\n")
        
        with open(output_file, 'a') as f:
            f.write(data_row + "\n")

    print('Completed in %s' % timedelta(seconds=(datetime.now() - start_time).seconds))
