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
import os
import sys

from log import info, warning

try:

    from arcpy import env, CheckOutExtension, AddField_management, CalculateField_management, MakeQueryTable_management, \
        CopyFeatures_management, FeatureClassToFeatureClass_conversion, CopyRaster_management
    from arcpy.sa import TopoPointElevation, TopoCliff, TopoToRaster, Extent, ExtractMultiValuesToPoints
except ImportError:
    warning('You are not using a version of python with arcpy installed')
    sys.exit(1)

def topo_to_raster(cfg, shapefile, i, elev_ID, unit):
    info('topo_to_raster', shapefile, i, elev_ID, unit)
    # ps_units = ['udockum_base_pecos', 'ldockum_base_pecos', 'uochoan_base_pecos', 'lochoan_base_pecos']
    # db_units = ['udockum_base', 'ldockum_base', 'uochoan_base', 'lochoan_base']

    # Set local variables
    in_point_elevations = TopoPointElevation([[shapefile, '{}'.format(elev_ID)]])

    in_cliff = TopoCliff(['cbp_f', 'gm_f'])

    in_features = ([in_point_elevations, in_cliff])

    # Check out the ArcGIS Spatial Analyst extension license
    CheckOutExtension("Spatial")

    # Execute TopoToRaster
    # arcpy.env.extent = r'E:\BASEMAP\dbasin_spadtm_resamp_ft.tif'
    outTTR = TopoToRaster(in_features, "1000", Extent(485690, 3530089, 693190, 3627589), "20", "#", "#", "ENFORCE",
                          "SPOT",
                          "20", "#", "1", "0", "0", "200", '#', '#', 'ERROR_FILE{}.txt'.format(unit), '#', '#', '#',
                          'ERROR_PTS_{}'.format(unit))

    name = 't2r{}{}.tif'.format(unit, i)
    p = os.path.join(cfg['output_directory'], name)
    outTTR.save(p)
    info(i, 't2r done')
    return p


def extract_vals_to_pts(shapefile, raster, i):
    # Name: ExtractMultiValuesToPoints_Ex_02.py
    # Description: Extracts the cells of multiple rasters as attributes in
    #    an output point feature class.  This example takes a multiband IMG
    #    and two GRID files as input.
    # Requirements: Spatial Analyst Extension

    # Set local variables
    in_point_features = shapefile
    in_raster_list = [[raster, 'ei{}'.format(i)]]

    # Check out the ArcGIS Spatial Analyst extension license
    CheckOutExtension("Spatial")

    # Execute ExtractValuesToPoints
    ExtractMultiValuesToPoints(in_point_features, in_raster_list)

    info(i, 'extract multi values to pts done')


def calc_sigma(shapefile, i, elev_id):
    name = f'sigmae_{i}'
    AddField_management(shapefile, name, "DOUBLE")

    name_abs = f'abs(!{elev_id}! - !ei{i}!)'
    CalculateField_management(shapefile, name,
                              name_abs, "PYTHON")
    info(i, 'calculate sigma done')


def query_by_sigma(shapefile, i, unit):
    table = shapefile
    out_table = f"{unit}_qc{i}_result_table"
    in_key_field_option = 'USE_KEY_FIELDS'
    in_key_field = ''
    in_field = ''
    where_clause = f'sigmae_{i} < 100'

    in_features = MakeQueryTable_management(table, out_table, in_key_field_option, in_key_field, in_field,
                                            where_clause)
    out_feature_class = f"{unit}_qc{i}_result"
    CopyFeatures_management(in_features, out_feature_class)

    print(i, 'query done and shapefile output')


def copy_master_to_working(unit):
    FeatureClassToFeatureClass_conversion(f'{unit}_data_all_orig', env.workspace,
                                          f'{unit}_data_all_working')

    working_fc = f'{unit}_data_all_working'
    return working_fc


def copy_final(cfg, raster, unit):
    out_geodatabase = cfg['out_geodb']
    # name = '{}\{}'.format(out_geodatabase, unit)
    name = os.path.join(out_geodatabase, unit)
    CopyRaster_management(raster, name, '', '', '-3.402823e38',
                          'NONE', 'NONE', '32_BIT_FLOAT', '', '')


def configure_arcpy(cfg):
    env.workspace = cfg['arcpy_workspace']
    env.snapRaster = cfg['arcpy_snapRaster']
    env.overwriteOutput = True


def do_step1(cfg):
    configure_arcpy(cfg)

    unitnames = ['alvbase',
                 'udockum_base',
                 'sr_top',
                 'ldockum_base',
                 'deweylake_base',
                 'uochoan_base',
                 'lochoan_base',
                 'artesia_base']
    elevIDs = ['alv_b_elev',
               'ud_b_elev',
               'srtop_elev',
               'ld_b_elev',
               'rtop_elev',
               'saltop_elev',
               'lochoan_b_elev',
               'art_b_elev']

    for unit, elev_ID in zip(unitnames, elevIDs):
        working_fc = copy_master_to_working(unit)

        temp_fcs = [working_fc,
                    f'{unit}_qc0_result',
                    f'{unit}_qc1_result',
                    f'{unit}_qc2_result',
                    f'{unit}_qc3_result']

        for i, fc in enumerate(temp_fcs):
            print(i, fc)
            out_ras = topo_to_raster(cfg, fc, i, elev_ID, unit)

            extract_vals_to_pts(fc, out_ras, i)
            calc_sigma(fc, i, elev_ID)
            query_by_sigma(fc, i, unit)
            if i == 4:
                copy_final(cfg, out_ras, unit)
            info(f'{unit} QC process finished')
# ============= EOF =============================================
