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
import pprint
import arcpy

import yaml

from log import warning


def load_configuration():
    p = 'config.yaml'
    if not os.path.isfile(p):
        p = os.path.join(os.getcwd(), p)
        warning(f'No Configuration file found at {p}')
        return

    with open(p, 'r') as rfile:
        return yaml.load(rfile, yaml.SafeLoader)


def validate_configuration(cfg):
    for check in (
             #'input_directory',
             # 'working_geodb_name',
              #    'output_directory',
                  'basin',):
        if check in cfg:
            continue
        warning(f'{check} missing from configuration')
        return

    return True


def setup_configuration(cfg):
    basin = cfg['basin']

    # create major workspace
    root_path = os.path.join('C:\\', 'geomod3d', basin, f"{basin}_geodata")
    if not os.path.exists(root_path):
        os.makedirs(root_path)
    cfg['basin_workspace'] = root_path

    # create base data folder with basin name
    basedata = os.path.join(cfg['basin_workspace'], 'a_basedata')
    if not os.path.exists(basedata):
        os.makedirs(basedata)
    cfg['base_folder'] = basedata

    # create base data geodatabase
    basedata_gdb = os.path.join(cfg['base_folder'], 'base.gdb')
    if not os.path.exists(basedata_gdb):
        arcpy.CreateFileGDB_management(basedata, 'base.gdb')
    cfg['base_gdb'] = basedata_gdb

    # create working data folder
    working_folder = os.path.join(root_path, 'b_working')
    if not os.path.exists(working_folder):
        os.makedirs(working_folder)
    cfg['working_folder'] = working_folder

    # create working gdb
    working_gdb = os.path.join(working_folder, 'working_geodb.gdb')
    if not os.path.exists(working_gdb):
        arcpy.CreateFileGDB_management(working_folder, 'working_geodb.gdb')
    cfg['working_geodb'] = working_gdb

    # create folder for stepwise output folders
    result_folder = os.path.join(root_path, 'c_results')
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    cfg['result_folder'] = result_folder

    # create stepwise output folders
    qc_out = os.path.join(result_folder, 'b001_qcresults')
    if not os.path.exists(qc_out):
        os.makedirs(qc_out)
    cfg['b001'] = qc_out

    modelras_out_f = os.path.join(result_folder, 'b002_modelsurfs')
    if not os.path.exists(modelras_out_f):
        os.makedirs(modelras_out_f)
    cfg['modelras_out_f'] = modelras_out_f

    modelras_out_g = os.path.join(modelras_out_f, 'b002.gdb')
    if not os.path.exists(modelras_out_g):
        arcpy.CreateFileGDB_management(modelras_out_f, 'b002.gdb')
    cfg['modelras_out_g'] = modelras_out_g

    mc_surfs = os.path.join(result_folder, 'b003_mc')
    if not os.path.exists(mc_surfs):
        os.makedirs(mc_surfs)
    cfg['b003'] = mc_surfs

    uncert_out = os.path.join(result_folder, 'b004_uncert')
    if not os.path.exists(uncert_out):
        os.makedirs(uncert_out)
    cfg['b004'] = uncert_out


def report_configuration(cfg):
    msg = ' Configuration '
    print(f'{msg:*^80}')
    print(pprint.pformat(cfg, width=10))
    print('*'*80)


# ============= EOF =============================================
