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
# Import system modules
import arcpy
from arcpy import env
from arcpy.sa import *
import os

# Set environment settings - USER SPECIFIC
# arcpy.env.snapRaster = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\BASEMAP.gdb\z_snapraster'
# arcpy.env.overwriteOutput = True
# env.outputCoordinateSystem = arcpy.SpatialReference("NAD 1983 UTM Zone 13N")
# arcpy.CheckOutExtension('Spatial')

# working_geodb = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\working_geodb_and_folders_ex\working_geodb_example.gdb'
# working_folder = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\working_geodb_and_folders_ex\working_folder_example'
# env.workspace = working_folder
# stepwise_output_folder = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\working_geodb_and_folders_ex\b001'

## USER INPUT - WORKSPACE - WHERE NETWORK GDB WILL BE CREATED (need to create this folder ahead of time)
# CONFIG_FIN_V = {'ROOT': r'W:\statewide\aqua_map3D\am_references\MarsMapMethods\tutorial\example_network_output'}

v = '001'


def init_arcpy(cfg):
    # arcpy.env.snapRaster = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\BASEMAP.gdb\z_snapraster'

    path = os.path.join(cfg['workingdata'], 'BASEMAP.gdb', 'z_snapraster')
    arcpy.env.snapRaster = path
    arcpy.env.overwriteOutput = cfg.get('overwriteOutput', True)
    env.outputCoordinateSystem = arcpy.SpatialReference("NAD 1983 UTM Zone 13N")
    arcpy.CheckOutExtension('Spatial')
    env.workspace = cfg['working_folder']


#
# def get_working_folder(cfg):
#     # working_folder = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\working_geodb_and_folders_ex\working_folder_example'
#     name = cfg.get('working_folder', 'working_folder_example')
#     return os.path.join(get_working_geodb_and_folders_ex(cfg), name)
#
#
# def get_workingdata(cfg):
#     basin = cfg['basin']
#     return os.path.join('E:', basin, f"{basin}_geodata", 'workingdata')
#
#
# # basinname = 'Delaware Basin'
# def get_working_geodb_and_folders_ex(cfg):
#     return os.path.join(get_workingdata(cfg), 'working_geodb_and_folders_ex')
#
#
# def get_working_geodb_path(cfg):
#     name = cfg.get('working_geodb_name', 'working_geodb_name_example')
#     if not name.endswith('.gdb'):
#         name = f'{name}.gdb'
#
#     working_geodb_and_folders_ex = get_working_geodb_and_folders_ex(cfg)
#     working_geodb = os.path.join(working_geodb_and_folders_ex, name)
#     return working_geodb


# def set_working_geodb_workspace(cfg):
#     working_geodb = get_working_geodb_path(cfg)
#     env.workspace = working_geodb

def get_extent(cfg):
    if cfg.get('extent'):
        cfg = cfg['extent']
        # ext = Extent(485690, 3530089, 693190, 3627589)
        ext = Extent(cfg['xmin'], cfg['ymin'], cfg['xmax'], cfg['ymax'])
    else:
        arcpy.env.extent = cfg['snapRaster']# r'E:\BASEMAP\dbasin_spadtm_resamp_ft.tif'
        ext = '#'
    return ext


def topo_to_raster(cfg, fc, i, elev_ID, unit):
    print(fc, i, elev_ID, unit)
    # set_working_geodb_workspace(cfg)
    env.workspace = cfg['working_geodb']

    # Set local variables
    inPointElevations = TopoPointElevation([[fc, '{}'.format(elev_ID)]])
    # inPointElevations = TopoPointElevation([[shapefile, '{}_elev'.format(unitname)]])
    # inBoundary = TopoBoundary([r'feature_classes\dbasin_bdry_buffer10km.shp'])
    # inContours = TopoContour([['contours.shp', 'spot_meter']])
    # inLake = TopoLake(['lakes.shp'])
    # inSinks = TopoSink([['sink1.shp', 'elevation'], ['sink2.shp', 'none']])
    # inStream = TopoStream(['streams.shp'])
    inCliff = TopoCliff(['cbp_f', 'gm_f'])
    # inCoast = TopoCoast(['coast.shp'])
    # inExclusion = TopoExclusion(['ignore.shp'])

    inFeatures = ([inPointElevations, inCliff])

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    # Execute TopoToRaster
    # arcpy.env.extent = r'E:\BASEMAP\dbasin_spadtm_resamp_ft.tif'

    ext = get_extent(cfg)

    outTTR = TopoToRaster(inFeatures, "1000", ext,
                          "20", "#", "#", "ENFORCE",
                          "SPOT",
                          "20", "#", "1", "0", "0", "200", '#', '#',
                          f'ERROR_FILE{unit}.txt',
                          '#', '#', '#',
                          f'ERROR_PTS_{unit}')

    name = f'b001002r{unit}{i}.tif'

    # stepwise_output_folder = os.path.join(get_working_geodb_and_folders_ex(cfg), 'b001')
    stepwise_output_folder = cfg['b001']
    path = os.path.join(stepwise_output_folder, name)
    outTTR.save(path)
    print(i, 't2r done')
    return path

    # b001002r = r'{}\b001002r{}{}.tif'.format(stepwise_output_folder, unit, i)
    # outTTR.save(b001002r)
    # print(i, 't2r done')
    # return b001002r

    # sys.exit('CHECK IF T2R WAS DONE CORRECTLY WITH CLIFF')
    # sys.exit('CHECK DIAGNOSTIC FILE')


def extract_vals_to_pts(cfg, fc, raster, i, unit):
    print('extracting values from', raster, 'to', fc)
    # Name: ExtractMultiValuesToPoints_Ex02.py
    # Description: Extracts the cells of multiple rasters as attributes in
    #    an output point feature class.  This example takes a multiband IMG
    #    and two GRID files as input.
    # Requirements: Spatial Analyst Extension
    env.workspace = cfg['working_geodb']
    # set_working_geodb_workspace(cfg)
    # Set local variables
    inPointFeatures = fc
    inRasterList = [[raster, f'ei{i}']]

    # Check out the ArcGIS Spatial Analyst extension license
    arcpy.CheckOutExtension("Spatial")

    # Execute ExtractValuesToPoints
    ExtractMultiValuesToPoints(inPointFeatures, inRasterList, 'NONE')

    print(i, 'extract multi values to pts done')


def calc_sigma(fc, i, elev_ID):
    arcpy.AddField_management(fc, f"sigmae_{i}", "DOUBLE")
    arcpy.CalculateField_management(fc, f"sigmae_{i}",
                                    f'abs(!{elev_ID}! - !ei{i}!)', "PYTHON")
    print(i, 'calculate sigma done')


def query_by_sigma(cfg, fc, i, unit):
    print('input file name = ', fc)
    process_tag = 'b001003'
    table = fc
    out_table = f'{process_tag}fc_{unit}_{i}_result_table'
    in_key_field_option = 'USE_KEY_FIELDS'
    in_key_field = ''
    in_field = ''
    where_clause = 'sigmae_{} < 100'.format(i)

    in_features = arcpy.MakeQueryTable_management(table, out_table, in_key_field_option, in_key_field, in_field,
                                                  where_clause)
    # working_folder = get_working_folder(cfg)
    out_shapefile = os.path.join(cfg['working_folder'],
                                 f'{process_tag}{unit}{i}.shp')

    # out_feature_class = r'{}\b001003{}{}'.format(working_geodb, unit, i)
    out_feature_class = os.path.join(cfg['working_geodb'], f'{process_tag}{unit}{i}')
    arcpy.CopyFeatures_management(in_features, out_shapefile)
    arcpy.CopyFeatures_management(in_features, out_feature_class)

    print(i, 'query done and shapefile output')


def copy_master_to_working(cfg, unit):
    # in_workspace = get_working_geodb_path(cfg)
    in_workspace = cfg['working_geodb']
    path = os.path.join(in_workspace, f'{unit}_data_all_orig')
    arcpy.FeatureClassToFeatureClass_conversion(path,
                                                in_workspace,
                                                f'{unit}_data_all_working')

    working_fc = '{}_data_all_working'.format(unit)

    return working_fc


def copy_final(cfg, raster, unit):
    # name = r'{}\b001002r{}.tif'.format(working_folder, unit)
    name = f'b001002r{unit}.tif'
    # working_folder = get_working_folder(cfg)
    name = os.path.join(cfg['working_folder'], name)
    arcpy.CopyRaster_management(raster, name, 'TIFF', '', '-3.402823e38',
                                'NONE', 'NONE', '32_BIT_FLOAT', '', '')

    return name


def ebm_modelbound(cfg, rasters, n, unitnames):
    print('BEGIN masking rasters to model boundary')
    # USER INPUT - STUDY BOUNDARY RASTER
    # mask = r'E:\3D_spatial_general\3d mapping areas\Delaware_Basin.shp'
    mask = cfg['model_extent_polygon']
    print('Model boundary raster or polygon used for clipping = {}'.format(mask))
    inRasters = rasters
    outRasters = []

    for raster, unitname in zip(inRasters, unitnames):
        print('...masking', raster, 'with', mask, 'unitname ==== {}'.format(unitname))
        outExtractByMask = ExtractByMask(raster, mask)
        name = r'b00200{}r{}.tif'.format(n, unitname)
        outExtractByMask.save(name)
        outRasters.append(name)
        print('Masked {} with {} and output file = {}'.format(raster, mask, name))

    print('FINISHED masking b00200{}r extents with model ext bdry'.format(n))
    return outRasters


def extract_by_mask(inRasters, maskdata, unitnames):
    print('BEGIN masking full extent rasters to discrete extent')
    print('inRasters = ', inRasters)
    print('inRasters[1:] = ', inRasters[1:])
    print('inRasters[0] = ', inRasters[0])
    inMaskData = maskdata
    disc_exts = ['{}'.format(inRasters[0])]
    inRasters = inRasters[1:]

    print('disc exts beginning = ', disc_exts)
    print('unitnames = ', unitnames)

    for raster, mask, unitname in zip(inRasters, inMaskData, unitnames[1:]):
        print('...masking', raster, 'with', mask, '= {}'.format(unitname))
        name = r'b002005r{}.tif'.format(unitname)
        outExtractByMask = ExtractByMask(raster, mask)
        outExtractByMask.save(name)
        disc_exts.append(name)
        print('disc_exts = ', disc_exts, 'after loop')
        print('...masked {} with {} ----- unitname = {}'.format(raster, mask, unitname))

    print('FINISHED masking full extent rasters to discrete extent')
    return disc_exts


def reclassify(rasters, unitnames):
    print('BEGIN reclassifying negative thickness values to NoData')
    reclassField = 'VALUE'
    inRasters = rasters
    remap = RemapRange([[-100000, 0, 'NoData']])

    reclass_rasters = []

    for ras, unitname in zip(inRasters, unitnames[1:]):
        print('...reclassifying', ras, '--- should end in --->', unitname)
        outReclass = Reclassify(ras, reclassField, remap)
        name = 'b002004r{}.tif'.format(unitname)
        outReclass.save(name)
        reclass_rasters.append(name)
        print('Reclassified negative thickness values for {} = {}'.format(ras, unitname))

    maskdata = reclass_rasters
    print('FINISHED reclassifying negative thickness values to NoData')
    return maskdata


def minus(rasters, unitnames):
    print('BEGIN calculating full extent unit thicknesses')
    # set_default_workspace()

    output_rasters = []

    for upper_unit, lower_unit, unitname in zip(rasters, rasters[1:], unitnames[1:]):
        print('...calculating', upper_unit, 'minus', lower_unit)
        outMinus = Minus(upper_unit, lower_unit)
        name = f'b002003r{unitname}.tif'
        outMinus.save(name)
        output_rasters.append(name)
        print('...calculated thickness of', lower_unit, '...should be same as...', unitname)

    thickness_ras = output_rasters
    print('FINISHED calculating full extent unit thicknesses')
    return thickness_ras


def mosaic_min(cfg, rasters, unitnames):
    print('BEGIN mosaicking upper/lower surfaces and taking minimum')
    fullext_rasters = [rasters[0]]

    # create alluvium base full extent first:
    name_alv_fe = f'b002002r{unitnames[1]}.tif'
    print('...mosaicking', rasters[0], 'with', rasters[1], 'MINIMUM', 'name', name_alv_fe)

    working_folder = cfg['working_folder']
    arcpy.MosaicToNewRaster_management(f'{rasters[0]};{rasters[1]}', working_folder,
                                       name_alv_fe, env.outputCoordinateSystem,
                                       '32_BIT_FLOAT', '1000', '1', 'MINIMUM', '')

    fullext_rasters.append(name_alv_fe)

    i = 1
    for lower_ras, unitname in zip(rasters[2:], unitnames[2:]):
        # now mosaic output fullext_raster with underlying unit
        fe_upper_ras = fullext_rasters[i]

        name_fullext = f'b002002r{unitname}.tif'
        print('...mosaicking', fe_upper_ras, 'with', lower_ras, 'MINIMUM', 'out_name', name_fullext)

        arcpy.MosaicToNewRaster_management(f'{fe_upper_ras};{lower_ras}', working_folder,
                                           name_fullext, env.outputCoordinateSystem,
                                           '32_BIT_FLOAT', '1000', '1', 'MINIMUM', '')
        fullext_rasters.append(name_fullext)
        i = i + 1

    mosaic_ras = fullext_rasters
    print('FINISHED mosaicking upper/lower surfaces and taking minimum')
    return mosaic_ras


def resample(surfaces, units, size='1000', kind='NEAREST'):
    # USER INPUT CELL SIZE
    size_label = '{}x{}m'.format(size, size)
    print('BEGIN resampling rasters to {}'.format(size_label))

    # resampled_rasters = ['{}_resample'.format(r) for r in surfaces]
    resampled_rasters = ['b002001r{}.tif'.format(u) for u in units]
    for surface, output in zip(surfaces, resampled_rasters):
        print('...resampling', surface)
        print('out name = ', output)
        arcpy.Resample_management(surface, output, size, kind)

    res_ras = resampled_rasters
    print('FINISHED resampling rasters to {}'.format(size_label))
    return res_ras


def place_final002_surfaces(cfg, rasters):
    outfolder002 = cfg['final_out_f']

    final_out_g = os.path.join(outfolder002, 'b002.gdb')
    if not os.path.exists(final_out_g):
        arcpy.management.CreateFileGDB(outfolder002, 'b002.gdb')
    outgdb002 = final_out_g

    for raster in rasters:
        print('raster = ', raster)
        # name = r'{}\{}'.format(outfolder002, raster)
        name = os.path.join(outfolder002, raster)
        print('output TIFF file = ', name)
        arcpy.CopyRaster_management(raster, name, '', '', '-3.402823e38', 'NONE', 'NONE',
                                    '32_BIT_FLOAT', '', '', 'TIFF', '')

        gdbname = os.path.join(outgdb002, raster.replace('.tif', ''))
        print('output GRID file = ', gdbname)
        arcpy.CopyRaster_management(raster, gdbname, '', '', '-3.402823e38', 'NONE', 'NONE',
                                    '32_BIT_FLOAT', '', '', 'GRID', '')


def copy_to_network_folder(cfg, versionrasters):
    print('COPYING VERSIONS TO NETWORK VERSION FOLDER')

    # output_location = r'{}'.format(CONFIG_FIN_V['ROOT'])
    output_location = cfg['network_output_directory']
    if not os.path.exists(output_location):
        os.makedirs(output_location)

    for raster in versionrasters:
        print('...copying', raster, 'to', output_location, '---- as ----', raster)
        # name = r'{}\v{}{}'.format(output_location, v, raster)
        name = os.path.join(output_location, f'v{v}{raster}')
        arcpy.CopyRaster_management(raster, name, '', '', '-3.402823e38',
                                    'NONE', 'NONE', '32_BIT_FLOAT', '', '')

    print('finished copying to network folder: {}'.format(output_location))


def make_network_gdb_name(cfg, name, create=True):
    basename = f'{name}.gdb'
    output_directory = cfg['network_output_directory']
    name = os.path.join(output_directory, basename)
    if create:
        # if does not exist create it
        if not arcpy.Exists(name):
            arcpy.CreateFileGDB_management(output_directory, basename)
    return name


def export_fc_to_csv(cfg, fc, elevID, unit):
    print(f'EXPORTING attribute table to csv for {fc}')
    fields = ['OBJECTID', 'Field1', 'OriginalID', 'Easting', 'Northing', f'{elevID}',
              'DataSource', 'ei4']
    out_path = cfg['working_folder']
    out_csv = f'{unit}_modelbuild_inputdata.csv'
    arcpy.ExportXYv_stats('{}'.format(fc), fields, 'COMMA', os.path.join(out_path, out_csv), 'ADD_FIELD_NAMES')
    print(f'EXPORTED attribute table to csv for {unit}')

    qc_csv = os.path.join(out_path, out_csv)

    return qc_csv


def uncert_topo_to_raster(cfg, name, n):
    # need to make n number of surfaces, where n = number of random sample text files generated in sample_percent

    env.workspace = cfg['working_folder']
    print('current workspace =', env.workspace)

    arcpy.CheckOutExtension("Spatial")

    for i in range(1, n + 1):
        print(f'Running topo to raster on {name}, {i}')
        # fc = f'{name}_fc00{i}.shp'
        fc = make_feature_class_name(name, i)
        inPointElevations = TopoPointElevation([[fc, 'ei4']])
        inCliffs = TopoCliff(['cbp_f.shp', 'gm_f.shp'])
        inFeatures = ([inPointElevations, inCliffs])

        ext = get_extent(cfg)

        outTTR = TopoToRaster(inFeatures, "1000",
                              #Extent(485690, 3530089, 693190, 3627589),
                              ext,
                              "20", "#", "#", "ENFORCE",
                              "SPOT",
                              "20", "#", "1", "0", "0", "200")

        pname = f'{name}_t2r{i}.tif'
        path = os.path.join(env.workspace, pname)
        print('out t2r path/name == ', path)
        outTTR.save(pname)
        print('t2r', f'{name}, {i}')
    print('t2r done')


def make_feature_class_name(name, i):
    return f'{name}_fc{i:03n}.shp'


def make_rand_csv_name(name, i):
    return f'{name}_rand{i:03n}.csv'


def create_feature_class(cfg, name, n):
    working_folder = cfg['working_folder']
    print('fc current workspace == ', env.workspace)
    print('working folder', working_folder)
    for i in range(1, n + 1):
        print(f'Creating feature class out of rand sample files: {name}, {i}')
        table = make_rand_csv_name(name, i)
        X = 'Easting'
        Y = 'Northing'
        Z = 'ei4'
        out_layername = f'{name}_lyr{i:03n}'
        out_layer = os.path.join(working_folder, out_layername)
        arcpy.MakeXYEventLayer_management(table, X, Y, out_layer, env.outputCoordinateSystem, Z)

        fc_name = make_feature_class_name(name, i)
        arcpy.FeatureClassToFeatureClass_conversion(out_layer, working_folder, fc_name)
        print('fc', f'{name} {i}', 'created')
    print('{} all feature classes created'.format(name))


def sample_percent(cfg, df, unitname, n):
    # ei4 is the final point dataset used in making the final raster surface, so need to take a random 75% of
    # that dataset (500) times - store in 'output_random_sample_datasets'
    # out_path = r'C:\Users\mfichera\PycharmProjects\3D_mapping\db_methods\randsample_exports'
    # out_path = working_folder
    # out_path = get_working_folder(cfg)
    out_path = cfg['working_folder']
    print('sample percent current workspace == ', env.workspace)
    print('outoyasdf ', out_path)
    for i in range(1, n + 1):
        print(f'RANDOM SAMPLE {i}: taking random sample of {unitname} dataset')
        rand_samp = df.sample(frac=.75)
        # p = f'{unitname}_rand00{i}.csv'
        p = make_rand_csv_name(unitname, i)
        rand_samp.to_csv(os.path.join(out_path, p))
        print(f'random sample file taken {i}, {unitname}')
    print('random sample files created')

    return out_path


# def find_qc_csvs():
#     qc_out_path = r'C:\Users\mfichera\PycharmProjects\3D_mapping\db_methods\qcexports'
#     qc_location = cfg['qc_location']
#     qc_name = cfg['qc_result_name']
#     qc_out_csv = f'{unit}_modelbuild_inputdata.csv'


def generate_uncertainty_maps(cfg, unitnames, n, masks):
    # set_workspace(cfg)
    # env.workspace = working_folder
    # env.workspace = get_working_folder(cfg)
    env.workspace = cfg['working_folder']
    umaps = []
    for unit, mask in zip(unitnames[1:], masks[1:]):
        print(f'generating uncertainty map for unit {unit}')
        rasterList = arcpy.ListRasters(f'{unit}_t2r*', 'TIF')
        print('==== raster list =====')
        print(rasterList)
        print('...calculating standard deviation at each cell...')
        outSTD = CellStatistics(rasterList, 'STD', 'DATA')
        print('...calculating mean at each cell...')
        outMEAN = CellStatistics(rasterList, 'MEAN', 'DATA')

        print('...clipping to geologic extent...')
        outSTD_map = ExtractByMask(outSTD, mask)
        outMEAN_map = ExtractByMask(outMEAN, mask)

        print('...saving result...')
        nameSTD = f'{unit}_stdev{n}.tif'
        outSTD_map.save(nameSTD)
        nameMEAN = f'{unit}_mean{n}.tif'
        outMEAN_map.save(nameMEAN)

        print('...saving results to separate folder...')

        # uncertainty_path = r'E:\DelawareBasin\DelawareBasin_geodata\workingdata\working_geodb_and_folders_ex\b004'
        # uncertainty_path = os.path.join(get_working_geodb_and_folders_ex(cfg), 'b004')
        uncertainty_path = cfg['b004']
        if not os.path.exists(uncertainty_path):
            os.makedirs(uncertainty_path)
        outSTD_map.save(os.path.join(uncertainty_path, nameSTD))
        outMEAN_map.save(os.path.join(uncertainty_path, nameMEAN))

        umaps.append(nameSTD)
        umaps.append(nameMEAN)

    return umaps
# ============= EOF =============================================
