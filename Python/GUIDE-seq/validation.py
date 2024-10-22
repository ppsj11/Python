"""
validation.py
=============

Contains utils for validating the filetype and existence of manifest-defined files/folders

"""

import logging
import os
import sys
from distutils.spawn import find_executable

logger = logging.getLogger('root')


def exists(filepath): # 파일이 존재하는지 확인.
    if not os.path.isfile(filepath): # 이 함수를 사용하여 주어진 경로에 파일이 존재하는지 체크. 파일이 존재하지 않으면 에러를 기록하고 프로그램 종료함.
        logger.error('{0} does not exist'.format(filepath))
        sys.exit()


def checkIfBinary(filepath): # 주어진 경로의 파일이 유효한 바이너리 실행파일인지 검사.
    executable = find_executable(filepath) # fine_executable로 실행파일을 찾고, 

    if executable is None: # 존재하지 않으면 에러를 기록.
        logger.error('Executable binary not found at {0}'.format(filepath))
        sys.exit()

    # First check if file exists
    exists(executable) # 파일이 존재하면,

    # Check if file is a valid binary
    # Adapted from http://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
    textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars)) # 파일의 처음 1024바이트를 읽어 유효한 바이너리인지 검사.

    if not is_binary_string(open(executable, 'rb').read(1024)):
        logger.error('{0} is not a valid binary'.format(executable))
        sys.exit()


def checkIfFasta(filepath): # FASTA 파일의 존재 유무 체크.
    # First check if file exists
    exists(os.path.abspath(filepath)) # 절대 결로를 이용하여 파일의 존재 체크


def checkIfFolder(folderpath): # 주어진 경로가 유효한 폴더인지 확인.
    # Check if the folder exists
    if not os.path.isdir(os.path.abspath(folderpath)): # os.path.isdir로 존재 여부를 체크하고, 존재하지 않으면 에러 기록.
        logger.error('{0} is not a valid folder path'.format(folderpath))
        sys.exit()


def checkIfValidUndemultiplexed(undemultiplexed): # 필요한 필드가 존재하는지 확인.
    # Check if read1, read2, index1, and index2 exist
    fields = ['forward', 'reverse', 'index1', 'index2']

    if set(fields) != set(undemultiplexed.keys()): # 특정 필드가 모든 존재하는지 체크하고 없으면 에러 기록.
        logger.error('Undemultiplexed field must contain references to "forward", "reverse", "index1", "index2"')
        sys.exit()

    invalid_file = False # 각 필드의 파일 존재 여부 체크. 유효하지 않은 파일이 있으면 에러 기록 후 종료.
    for field in fields:
        if not os.path.isfile(undemultiplexed[field]):
            logger.error('"read1" undemultiplexed field does not reference a valid file')
            invalid_file = True

    if invalid_file:
        sys.exit()


def checkIfValidSamples(samples): # 샘플 데이터가 유효한지 검증.
    # Check if control is one of the samples
    if 'Control' not in samples: # Control 샘플이 있는지 확인하고, 샘플이 정의되지 않은 경우 에러 기록.
        logger.error('A Control sample must be specified')
        sys.exit()

    if len(samples.keys()) == 0:
        logger.error('No samples defined')
        sys.exit()

    for sample in samples: # 각 샘플에 대해 필요한 필드가 모두 존재하는지 확인.
        if 'barcode1' not in samples[sample] or 'barcode2' not in samples[sample]:
            logger.error('barcode1 and barcode2 must be specified for {0} sample'.format(sample))
            sys.exit()
        if 'target' not in samples[sample]:
            logger.error('target sequence must be specified for {0} sample'.format(sample))
            sys.exit()


def validateManifest(manifest_data): # manifest에 필수 필드가 모두 존재하는지 검증.
    # Check if manifest contains the required fields
    fields = ['bwa', 'bedtools', 'reference_genome', 'output_folder', 'samples', 'undemultiplexed']
    missing_fields = False

    for field in fields:
        if field not in manifest_data.keys():
            logger.error('"{0}" field must be specified in manifest'.format(field))
            missing_fields = True

    if missing_fields:
        sys.exit()

    # Now validate each field
    checkIfBinary(manifest_data['bwa'])
    checkIfBinary(manifest_data['bedtools'])
    checkIfFasta(manifest_data['reference_genome'])
    checkIfValidUndemultiplexed(manifest_data['undemultiplexed'])
    checkIfValidSamples(manifest_data['samples'])
