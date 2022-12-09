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
import contextlib


@contextlib
def displayblock(txt, pad='~'):
    try:
        tag = f'{txt} started'
        paddedmessage(tag, pad)
        yield
    finally:
        tag = f'{txt} finished'
        paddedmessage(tag, pad)


def paddedmessage(tag, pad):
    print('~' * 80)
    if '|' in tag:
        tag = tag.split('|')
    if not isinstance(tag, (tuple, list)):
        tag = (tag,)

    for ti in tag:
        ti = f' {ti} '
        print(f'{ti:{pad}^80}')
    print('~' * 80)


def warning(msg):
    print(f'WARNING: ************ {msg}')

# ============= EOF =============================================
