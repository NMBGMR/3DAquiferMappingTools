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
             'working_geodb_name',
              #    'output_directory',
                  'basin',):
        if check in cfg:
            continue
        warning(f'{check} missing from configuration')
        return

    return True


def setup_configuration(cfg):
    basin = cfg['basin']
    workingdata = os.path.join('E:', basin, f"{basin}_geodata", 'workingdata')
    cfg['workingdata'] = workingdata

    working_geodb_and_folders = os.path.join(workingdata, 'working_geodb_and_folders_ex')
    cfg['working_geodb_and_folders'] = working_geodb_and_folders

    working_folder = os.path.join(working_geodb_and_folders, cfg['working_folder'])
    cfg['working_folder'] = working_folder

    working_geodb = os.path.join(working_folder, cfg['working_geodb_name'])
    cfg['working_geodb'] = working_geodb

    cfg['b001'] = os.path.join(working_geodb_and_folders, 'b001')
    cfg['b004'] = os.path.join(working_geodb_and_folders, 'b004')


def report_configuration(cfg):
    msg = ' Configuration '
    print(f'{msg:*^80}')
    print(pprint.pformat(cfg, width=10))
    print('*'*80)


# ============= EOF =============================================
