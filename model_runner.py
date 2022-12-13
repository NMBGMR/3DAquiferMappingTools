# ===============================================================================
# Copyright 2022 ross
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ===============================================================================

## written by marissa m. fichera

# Import system modules

import arcpy
from arcpy import env
from arcpy.sa import *
import os
import pandas as pd
import sys

from config import load_configuration, report_configuration, validate_configuration, setup_configuration
from log import paddedmessage, displayblock,welcome
from model import copy_master_to_working, topo_to_raster, extract_vals_to_pts, copy_final, export_fc_to_csv, calc_sigma, \
    query_by_sigma, resample, mosaic_min, minus, reclassify, extract_by_mask, ebm_modelbound, place_final002_surfaces, \
    sample_percent, create_feature_class, uncert_topo_to_raster, generate_uncertainty_maps, copy_to_network_folder, \
    init_arcpy, copy_basedata


def main():
    welcome()
    cfg = load_configuration()

    if not validate_configuration(cfg):
        report_configuration(cfg)
        return

    setup_configuration(cfg)
    report_configuration(cfg)


    init_arcpy(cfg)



    basinname = cfg['basin']
    unitnames = cfg['unitnames']
    faults = cfg['fault_list']
    elevIDs = cfg['elevIDs']

    copy_basedata(cfg)



    qc_output_rasters = []
    qc_output_fcs = []
    qc_output_csvs = []
    if cfg.get('do_qc', True):
        for unit, elev_ID in zip(unitnames, elevIDs):
            b001001fc = copy_master_to_working(cfg, unit)

            temp_fcs = [b001001fc,
                        f'b001003{unit}0',
                        f'b001003{unit}1',
                        f'b001003{unit}2',
                        f'b001003{unit}3']
            ntemp_fcs = len(temp_fcs)
            for i, fc in enumerate(temp_fcs):
                print(f'iteration {i}, {fc}')
                b001002r = topo_to_raster(cfg, fc, i, elev_ID, unit, faults)
                extract_vals_to_pts(cfg, fc, b001002r, i, unit)
                if i == (ntemp_fcs-1):
                    qc_out_r = copy_final(cfg, b001002r, unit)
                    qc_output_rasters.append(qc_out_r)
                    qc_output_fcs.append(fc)
                    qc_csv = export_fc_to_csv(cfg, fc, elev_ID, unit)
                    qc_output_csvs.append(qc_csv)
                    print('{} QC process finished'.format(unit))
                else:
                    calc_sigma(fc, i, elev_ID)
                    query_by_sigma(cfg, fc, i, unit)
    else:
        print('Skipping QC process')
        print('*********Make sure .csv files containing output QC data are in working folder*********')
        for unit in unitnames:
            out_csv = f'{unit}_modelbuild_inputdata.csv'
            qc_csv = os.path.join(cfg['working_folder'], out_csv)
            qc_output_csvs.append(qc_csv)

        print('DONE WITH QC FOR ALL UNITS')

    print(f'Raster surfaces going into model build = {qc_output_rasters}')

    paddedmessage('BUILDING|MODEL|STARTS|NOW')

    ## Step 02 - build down model
    b002_rasters = qc_output_rasters
    b002_rasters.insert(0, 'DEM.tif')
    unitnames.insert(0, 'DEM')

    env.workspace = cfg['working_folder']
    print('unitnames ====== ', unitnames)
    with displayblock('RESAMPLE RASTERS', pad='='):
        res_ras = resample(b002_rasters, unitnames)
        print(res_ras)

    with displayblock('FULL EXT. RASTERS'):
        fe_ras = mosaic_min(cfg, res_ras, unitnames)
        print(fe_ras)

    with displayblock('MINUS RASTERS'):
        thickness_ras = minus(fe_ras, unitnames)
        print(thickness_ras)

    with displayblock("MASK RASTERS"):
        maskdata = reclassify(thickness_ras, unitnames)
        print(maskdata)

    with displayblock("DISCRETE EXTENT RASTERS"):
        disc_exts = extract_by_mask(fe_ras, maskdata, unitnames)
        print(disc_exts)

    with displayblock('CLIPPING FULL AND DISCRETE EXTENT RASTERS TO MODEL EXTENT'):
        #### clips to model extent #####
        paddedmessage('FULL EXTENT')
        exts_modelext_full = ebm_modelbound(cfg, fe_ras[1:], '6', unitnames[1:])

        paddedmessage('DISCRETE EXTENT')
        # print('>>>>>>>>>>>>>> DISCRETE EXTENT <<<<<<<<<<<<<<<<<<')
        exts_modelext_disc = ebm_modelbound(cfg, disc_exts[1:], '7', unitnames[1:])

        paddedmessage('ISOPACH MAPS')
        isopach_modelext = ebm_modelbound(cfg, maskdata, '8', unitnames[1:])

    with displayblock('placing final discrete extent, full extent, and isopach rasters in separate local folder and '
                      'gdb', pad='+'):
        place_final002_surfaces(cfg, exts_modelext_disc)
        place_final002_surfaces(cfg, exts_modelext_full)
        place_final002_surfaces(cfg, isopach_modelext)

    with displayblock('generate uncertainty datasets'):
        n = 10
        if not qc_output_csvs:
            sys.exit('ERROR: QC CSV files are not in working folder')
        for c, unit in zip(qc_output_csvs, unitnames[1:]):
            print('output_csv = ', c)
            print('unitname = ', unit)
            qcdata = pd.read_csv(c)
            df = pd.DataFrame(qcdata)

            sample_percent(cfg, df, unit, n)
            create_feature_class(cfg, unit, n)
            uncert_topo_to_raster(cfg, unit, n)

    with displayblock('generate uncertainty maps'):
        umaps = generate_uncertainty_maps(cfg, unitnames, n, exts_modelext_disc)

    ########## copies to network ##############
    ########## NMBG specific and user input req'd ###########

    # type 'yes' as copy_command between the quotes if you want to copy to network
    # type 'no' between the quotes if you don't
    # copy_command = 'no'
    copy_command = input('Do you want me to copy the final surfaces to the network? [y]/n: ')
    # print('Do you want me to copy the final surfaces to the network? user says: ', copy_command)
    copy_command = copy_command.lower()
    if not copy_command or copy_command == 'y':
        print('okay - COPYING TO NETWORK FOLDER')
        copy_to_network_folder(cfg, exts_modelext_disc)
        copy_to_network_folder(cfg, exts_modelext_full)
        copy_to_network_folder(cfg, isopach_modelext)
        copy_to_network_folder(cfg, umaps)
        print('COPIED TO NETWORK FOLDER')
    elif copy_command == 'n':
        print('okay - DID NOT COPY TO NETWORK')
    else:
        print('error: please specify whether or not to copy to a network folder with y or n')

    paddedmessage(f'{basinname}|SUPERMODEL|IS|NOW|COMPLETE')


if __name__ == '__main__':
    main()

# ============= EOF =============================================
