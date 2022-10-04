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
from datetime import datetime

import click
import yaml


def warning(txt):
    click.secho(txt, fg='red')


def load_configuration():
    p = 'config.yaml'
    if not os.path.isfile(p):
        p = os.path.join(os.getcwd(), p)
        warning(f'No Configuration file found at {p}')
        return

    with open(p, 'r') as rfile:
        return yaml.load(rfile, yaml.SafeLoader)


def validate_configuration(cfg):
    for check in ('input_directory', 'output_directory'):
        if check in cfg:
            continue
        warning(f'{check} missing from configuration')
        return

    return True


def report_configuration(cfg):
    click.secho(pprint.pformat(cfg, width=10), fg='yellow')


def make_output_directory(cfg):
    out = cfg['output_directory']
    if not cfg.get('overwrite_output') and os.path.isdir(out):
        dt = datetime.now().strftime('%Y%m%d_%H%M%S')
        out = f'{out}{dt}'

    if not os.path.isdir(out):
        os.mkdir(out)


def do_model():
    cfg = load_configuration()
    report_configuration(cfg)
    if not validate_configuration(cfg):
        return

    make_output_directory(cfg)


if __name__ == '__main__':
    do_model()
# ============= EOF =============================================
