#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

guideseq.py
===========
serves as the wrapper for all guideseq pipeline

"""

import os # operating system: 운영체제에서 제공되는 여러 기능 수행 (파일 복사, 디렉터리 생성, 파일 목록 등)
import sys # system-specific parameters and functions: 인터프리터를 제어하는 모듈
import yaml # YAML Ain't Markup Language : 사람이 읽을 수 있는 데이터 직렬화 언어
import argparse # 스크립트 실행시 사용자료부터 입력을 받거나, 다양한 옵션 설정시 사용
import traceback # 오류 발생 시 상세한 오류 정보를 출력

# Set up logger
import log # 프로그램의 실행 과정, 상태, 오류 등을 기록하기 위해 사용
logger = log.createCustomLogger('root') # root 로거 생성

from alignReads import alignReads # alignReads 모듈에서 alignReads라는 이름의 함수 또는 클래스를 가져옴
from filterBackgroundSites import filterBackgroundSites
from umi import demultiplex, umitag, consolidate
from visualization import visualizeOfftargets
import identifyOfftargetSites
import validation

DEFAULT_DEMULTIPLEX_MIN_READS = 10000
DEFAULT_WINDOW_SIZE = 25
DEFAULT_MAX_SCORE = 7

CONSOLIDATE_MIN_QUAL = 15 # 최소 품질 점수 : 15
CONSOLIDATE_MIN_FREQ = 0.9 # 최소 빈도 : 특정 염기 또는 서열이 전체 데이터에서 90% 이상이여야 함

# class : 사용자 정의 데이터 구조를 정의하는 데 사용
# 객체: 클래스의 인스턴스로, 클래스에서 정의한 속성(데이터)과 메서드(함수)를 가지는 구체적인 실체
# 속성: 클래스 내에서 정의된 변수로, 객체의 상태를 나타냄
# 메서드: 클래스 내에서 정의된 함수로, 객체의 행동을 정의

class GuideSeq:

    def __init__(self):
        pass

    def parseManifest(self, manifest_path): # parseManifest라는 이름의 메서드 정의, self는 현재 인스턴스 참조, manifest_path는 매니페스트 파일의 경로를 나타내는 매개변수
        logger.info('Loading manifest...') # 로깅 시스템을 사용하여 "Loading manifest..."라는 정보를 로그로 기록. 매니페스트 파일을 로드하고 있다는 상태를 나타냄

        with open(manifest_path, 'r') as f: # with 구문을 사용하여 manifest_path에서 파일을 읽기 모드 'r'로 오픈. with 구문은 파일을 안전하게 열고, 작업이 끝난 후 자동으로 파일을 닫아주는 역할. f는 열린 파일 개체를 참조.
            manifest_data = yaml.safe_load(f) # yaml 모듈의 safe_load 함수를 사용하여 파일의 내용을 YAML 형식으로 파싱하고 그 결과는 manifest_data라는 변수에 저장됨. safe_load는 안전하게 YAML 데이터를 로드하는 방법으로, 외부 코드 실행을 방지함.
        
        if not "cores" in manifest_data: # manifest_data에 "cores"라는 키가 포함되어 있는지 확인. 만약 포함되어 있지 않다면, 다음 코드를 실행함.
            manifest_data['cores'] = 4 # 만약 "cores"키가 없다면, 기본값으로 4를 설정. 이는 manifest가 CPU 코어 수를 지정하지 않았을 경우, 기본값으로 4를 사용하독 하는 것임.
        
        # Set default tag/primer sequences if not specified
        if not "primer1" in manifest_data:
            manifest_data['primer1'] = 'TTGAGTTGTCATATGTTAAT'
        if not "primer2" in manifest_data:
            manifest_data['primer2'] = 'ACATATGACAACTCAATTAA'

        try: # try 블록: 예외 처리를 위해 사용. 일반적으로 프로그램 실행 중에 발생할 수 있는 오류를 사전에 처리하여, 프로그램이 비정상적으로 종료되는 것을 방지하고, 사용자에게 유용한 정보를 제공하기 위해 사용. validation.validateManifes(manifest_data) 함수를 호출.
            # Validate manifest data
            validation.validateManifest(manifest_data) # manifest_data의 유효성을 검증. 만약 유효성 검증에서 오류가 발생하면, 예외를 처리할 수 있는 구조.

            self.cores = manifest_data['cores'] # 유효성 검증이 성공하면, manifest_data에서 필요한 정보를 클래스의 속성으로 할당함. 각 속성은 manifest파일에서 읽어온 값으로 초기화됨.
            self.BWA_path = manifest_data['bwa']
            self.bedtools = manifest_data['bedtools']
            self.reference_genome = manifest_data['reference_genome']
            self.output_folder = manifest_data['output_folder']
            self.undemultiplexed = manifest_data['undemultiplexed']
            self.samples = manifest_data['samples']
            self.primer1 = manifest_data['primer1']
            self.primer2 = manifest_data['primer2']

        except Exception as e: # try 블록에서 발생한 모든 예외를 포착. Exception은 모든 일반적인 예외의 기본 클래스이므로, 이 구문은 어던 종류의 예외가 발생하더라고 처리할 수 있담. e는 발생한 예외 객체를 참조함.
            logger.error('Incorrect or malformed manifest file. Please ensure your manifest contains all required fields.') # logger 객체를 사용하여 오류 메시지를 기록. 
            sys.exit() # 프로그램을 종료하는 명령. 비정상 종료. 사용자에게 문제를 알리고 프로그램이 계속 진행되지 않도록 함.

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS
        # Allow the user to specify window size for off-target search
        if 'search_radius' in manifest_data:
            self.search_radius = manifest_data['search_radius']
        else:
            self.search_radius = DEFAULT_WINDOW_SIZE
        # Allow the user to specify window size for off-target search
        if 'max_score' in manifest_data:
            self.max_score = manifest_data['max_score']
        else:
            self.max_score = DEFAULT_MAX_SCORE
        # Allow the user to specify PAM seq. Yichao 3/6/2020
        if 'PAM' in manifest_data:
            self.PAM = manifest_data['PAM']
        else:
            self.PAM = "NGG"

        # Make sure the user has specified a control barcode
        if 'control' not in self.samples.keys():
            raise AssertionError('Your manifest must have a control sample specified.')

        # Make sure the user has both a sample and a control
        if len(self.samples) < 2:
            raise AssertionError('Your manifest must have at least one control and one treatment sample.')

        logger.info('Successfully loaded manifest.')

    def parseManifestDemultiplex(self, manifest_path): # parseManifestDemultiplex는 인스턴스 메서드로, minifest_path라는 인자를 받음.
        logger.info('Loading manifest for demultiplexing...') # 매니페스트 파일을 로드하고 있다는 정보를 로그로 기록.

        with open(manifest_path, 'r') as f: # with 구문을 사용하여 minifest_path에 있는 파일을 읽기 모드로 오픔. f는 열린 파일 객체.
            manifest_data = yaml.load(f) # yaml.load(f)를 사용하여 파일의 내용을 파싱함. manifest_data는 YAML 데이터가 저장된 사전임.

            try: # try 블록 내에서 매니페스트 데이터에서 필요한 필드를 가져옴.
                self.output_folder = manifest_data['output_folder'] # manifest_data에서 각 필드 값을 각각 클래스 인스턴스의 속성으로 할당함.
                self.undemultiplexed = manifest_data['undemultiplexed']
                self.samples = manifest_data['samples']

            except Exception as e:# 만약 위의 try 블록에서 오류가 발생하면, except 블록 실행됨.
                logger.error('Incomplete or incorrect manifest file. Please ensure your manifest contains all required fields.') # 오류 메시지를 로그로 기록하고, 사용자에게 알림.
                quit() # 프로그램을 종료

        # Allow the user to specify min reads for demultiplex if they want
        if 'demultiplex_min_reads' in manifest_data:
            self.demultiplex_min_reads = manifest_data['demultiplex_min_reads']
        else:
            self.demultiplex_min_reads = DEFAULT_DEMULTIPLEX_MIN_READS

        logger.info('Successfully loaded manifest for single-step demultiplexing.')

    def demultiplex(self):

        logger.info('Demultiplexing undemultiplexed files...')

        # Take our two barcodes and concatenate them
        swapped_sample_barcodes = {} # 사전 초기화. 사전에 생성된 바코드와 샘플 이름을 매핑하기 위해 사용됨.
        for sample in self.samples: # self.samples라는 샘플 데이터 구조를 반복함. 여기서 sample은 각 샘플의 키를 나타냄.
            barcode1 = self.samples[sample]['barcode1']
            barcode2 = self.samples[sample]['barcode2'] # 각 샘플의 바코드롤 가져옴.
            barcode = barcode1[1:8] + barcode2[1:8] # 두 barcode의 1-7번 인덱스까지 잘라낸 후 결합하여 새로운 바코드 생성. 
            swapped_sample_barcodes[barcode] = sample # 생성된 barcode를 키로 사용하고, 해당 샘플 이름을 값으로 하여 swapped_sample_barcodes 사전에 추가함.

        try: # try 블록 내에서 demultiplex 메서드 호출.
            demultiplex.demultiplex(self.undemultiplexed['forward'],
                                    self.undemultiplexed['reverse'],
                                    self.undemultiplexed['index1'],
                                    self.undemultiplexed['index2'],
                                    swapped_sample_barcodes,
                                    os.path.join(self.output_folder, 'demultiplexed'),
                                    min_reads=self.demultiplex_min_reads)

            self.demultiplexed = {}
            for sample in self.samples:
                self.demultiplexed[sample] = {}
                self.demultiplexed[sample]['read1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r1.fastq')
                self.demultiplexed[sample]['read2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.r2.fastq')
                self.demultiplexed[sample]['index1'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i1.fastq')
                self.demultiplexed[sample]['index2'] = os.path.join(self.output_folder, 'demultiplexed', sample + '.i2.fastq')

            logger.info('Successfully demultiplexed reads.')
        except Exception as e:
            logger.error('Error demultiplexing reads.')
            logger.error(traceback.format_exc())
            quit()

    def umitag(self):
        logger.info('umitagging reads...')

        try:
            self.umitagged = {}
            for sample in self.samples:
                self.umitagged[sample] = {}
                self.umitagged[sample]['read1'] = os.path.join(self.output_folder, 'umitagged', sample + '.r1.umitagged.fastq')
                self.umitagged[sample]['read2'] = os.path.join(self.output_folder, 'umitagged', sample + '.r2.umitagged.fastq')

                umitag.umitag(self.demultiplexed[sample]['read1'],
                              self.demultiplexed[sample]['read2'],
                              self.demultiplexed[sample]['index1'],
                              self.demultiplexed[sample]['index2'],
                              self.umitagged[sample]['read1'],
                              self.umitagged[sample]['read2'],
                              os.path.join(self.output_folder, 'umitagged'))

            logger.info('Successfully umitagged reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def consolidate(self, min_freq=CONSOLIDATE_MIN_FREQ, min_qual=CONSOLIDATE_MIN_QUAL):
        logger.info('Consolidating reads...')

        try:
            self.consolidated = {}

            for sample in self.samples:
                self.consolidated[sample] = {}
                self.consolidated[sample]['read1'] = os.path.join(self.output_folder, 'consolidated', sample + '.r1.consolidated.fastq')
                self.consolidated[sample]['read2'] = os.path.join(self.output_folder, 'consolidated', sample + '.r2.consolidated.fastq')

                consolidate.consolidate(self.umitagged[sample]['read1'], self.consolidated[sample]['read1'], min_qual, min_freq)
                consolidate.consolidate(self.umitagged[sample]['read2'], self.consolidated[sample]['read2'], min_qual, min_freq)

            logger.info('Successfully consolidated reads.')
        except Exception as e:
            logger.error('Error umitagging')
            logger.error(traceback.format_exc())
            quit()

    def alignReads(self):
        logger.info('Aligning reads...')

        try:
            self.aligned = {}
            for sample in self.samples:
                sample_alignment_path = os.path.join(self.output_folder, 'aligned', sample + '.sam')
                alignReads(self.cores,
                           self.BWA_path,
                           self.reference_genome,
                           self.consolidated[sample]['read1'],
                           self.consolidated[sample]['read2'],
                           sample_alignment_path)
                self.aligned[sample] = sample_alignment_path
                logger.info('Finished aligning reads to genome.')

        except Exception as e:
            logger.error('Error aligning')
            logger.error(traceback.format_exc())
            quit()

    def identifyOfftargetSites(self):
        logger.info('Identifying offtarget sites...')

        try:
            self.identified = {}

            # Identify offtarget sites for each sample
            for sample in self.samples:

                # Prepare sample annotations
                sample_data = self.samples[sample]
                annotations = {}
                annotations['Description'] = sample_data['description']
                annotations['Targetsite'] = sample

                if sample == 'control':
                    annotations['Sequence'] = ''
                else:
                    annotations['Sequence'] = sample_data['target']

                samfile = os.path.join(self.output_folder, 'aligned', sample + '.sam')

                self.identified[sample] = os.path.join(self.output_folder, 'identified', sample + '_identifiedOfftargets.txt')

                identifyOfftargetSites.analyze(samfile, self.reference_genome, self.identified[sample], annotations,
                                               self.search_radius, self.max_score, self.primer1, self.primer2)

            logger.info('Finished identifying offtarget sites.')

        except Exception as e:
            logger.error('Error identifying offtarget sites.')
            logger.error(traceback.format_exc())
            quit()

    def filterBackgroundSites(self):
        logger.info('Filtering background sites')

        try:
            self.filtered = {}

            # Filter background in each sample
            for sample in self.samples:
                if sample != 'control':
                    self.filtered[sample] = os.path.join(self.output_folder, 'filtered', sample + '_backgroundFiltered.txt')
                    filterBackgroundSites(self.bedtools, self.identified[sample], self.identified['control'], self.filtered[sample])
                    logger.info('Finished background filtering for {0} sample'.format(sample))

            logger.info('Finished filtering background sites.')

        except Exception as e:
            logger.error('Error filtering background sites.')
            logger.error(traceback.format_exc())

    def visualize(self):
        logger.info('Visualizing off-target sites')

        # try:
            # for sample in self.samples:
                # if sample != 'control':
                    # infile = self.identified[sample]
                    # outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')
                    # visualizeOfftargets(infile, outfile, title=sample)

            # logger.info('Finished visualizing off-target sites')

        # except Exception as e:
            # logger.error('Error visualizing off-target sites.')
            # logger.error(traceback.format_exc())

        for sample in self.samples: ## 3/6/2020 Yichao solved: visualization stopped when one sample failed
            if sample != 'control':
                try:
                    infile = self.identified[sample]
                    outfile = os.path.join(self.output_folder, 'visualization', sample + '_offtargets')
                    try:
                        self.PAM
                        visualizeOfftargets(infile, outfile, title=sample,PAM=self.PAM)
                    except:
                        visualizeOfftargets(infile, outfile, title=sample,PAM="NGG")
                except Exception as e:
                    logger.error('Error visualizing off-target sites: %s'%(sample))
                    logger.error(traceback.format_exc())
        logger.info('Finished visualizing off-target sites')

def parse_args(): # parse_args라는 함수 정의
    parser = argparse.ArgumentParser() # ArgumentParser 객체 생성 argparse.ArgumentParser()를 사용하여 인수 파서를 생성.

    subparsers = parser.add_subparsers(description='Individual Step Commands', # 하위 명령어 추가. add_subparsers 메서드를 사용하여 여러 하위 명령어를 추가할 수 있는 객체를 생성.
                                       help='Use this to run individual steps of the pipeline', # 각 하위 명령어에 대한 설명과 도움말을 설정함.
                                       dest='command') # 하위 명령어의 이름이 args.command로 저장됨.
    
    # 각 하위 명령어에 대해 필요한 인수와 도움말을 추가함.
    all_parser = subparsers.add_parser('all', help='Run all steps of the pipeline') # 모든 단계의 파이프라인을 실햄함.
    all_parser.add_argument('--manifest', '-m', help='Specify the manifest Path', required=True) # 매니페스트 파일의 경로를 요구하는 필수 인수
    all_parser.add_argument('--identifyAndFilter', action='store_true', default=False) # 이 플래그가 설정되면 식별 및 필터링 단계를 수행하도록 지정
    all_parser.add_argument('--skip_demultiplex', action='store_true', default=False) # 이 플래그가 설정되면 demultiplexing 단계를 건너뜀.

    demultiplex_parser = subparsers.add_parser('demultiplex', help='Demultiplex undemultiplexed FASTQ files') # fastq파일을 demultiplexing 함.
    demultiplex_parser.add_argument('--manifest', '-m', help='Specify the manifest path', required=True) # 매니페스트 파일의 경로를 요구하는 필수 인수

    umitag_parser = subparsers.add_parser('umitag', help='UMI tag demultiplexed FASTQ files for consolidation')
    umitag_parser.add_argument('--read1', required=True)
    umitag_parser.add_argument('--read2', required=True)
    umitag_parser.add_argument('--index1', required=True)
    umitag_parser.add_argument('--index2', required=True)
    umitag_parser.add_argument('--outfolder', required=True)

    consolidate_parser = subparsers.add_parser('consolidate', help='Consolidate UMI tagged FASTQs')
    consolidate_parser.add_argument('--read1', required=True)
    consolidate_parser.add_argument('--read2', required=True)
    consolidate_parser.add_argument('--outfolder', required=True)
    consolidate_parser.add_argument('--min_quality', required=False, type=float)
    consolidate_parser.add_argument('--min_frequency', required=False, type=float)

    align_parser = subparsers.add_parser('align', help='Paired end read mapping to genome')
    align_parser.add_argument('--bwa', required=True)
    align_parser.add_argument('--genome', required=True)
    align_parser.add_argument('--read1', required=True)
    align_parser.add_argument('--read2', required=True)
    align_parser.add_argument('--outfolder', required=True)

    identify_parser = subparsers.add_parser('identify', help='Identify GUIDE-seq offtargets')
    identify_parser.add_argument('--aligned', required=True)
    identify_parser.add_argument('--genome', required=True)
    identify_parser.add_argument('--outfolder', required=True)
    identify_parser.add_argument('--target_sequence', required=True)
    identify_parser.add_argument('--description', required=False)
    identify_parser.add_argument('--max_score', required=False, type=int, default=7)
    identify_parser.add_argument('--search_radius', required=False, type=int, default=25)

    filter_parser = subparsers.add_parser('filter', help='Filter identified sites from control sites')
    filter_parser.add_argument('--bedtools', required=True)
    filter_parser.add_argument('--identified', required=True)
    filter_parser.add_argument('--background', required=True)
    filter_parser.add_argument('--outfolder', required=True)

    visualize_parser = subparsers.add_parser('visualize', help='Visualize off-target sites')
    visualize_parser.add_argument('--infile', required=True)
    visualize_parser.add_argument('--outfolder', required=True)
    visualize_parser.add_argument('--title', required=False)

    return parser.parse_args()


def main(): # main이라는 함수 정의. 프로그램의 주요 실행 흐름을 담당.
    args = parse_args()

    if args.command == 'all':

        if args.identifyAndFilter:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)

                # Bootstrap the aligned samfile paths
                g.aligned = {}
                for sample in g.samples:
                    g.aligned[sample] = os.path.join(g.output_folder, 'aligned', sample + '.sam')

                g.identifyOfftargetSites()
                g.filterBackgroundSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        elif args.skip_demultiplex:
            try:
                g = GuideSeq()
                g.parseManifest(args.manifest)
                g.demultiplexed = {}
                for sample in g.samples:
                    g.demultiplexed[sample] = {}
                    g.demultiplexed[sample]['read1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r1.fastq')
                    g.demultiplexed[sample]['read2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.r2.fastq')
                    g.demultiplexed[sample]['index1'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i1.fastq')
                    g.demultiplexed[sample]['index2'] = os.path.join(g.output_folder, 'demultiplexed', sample + '.i2.fastq')
                    if not os.path.isfile(g.demultiplexed[sample]['read1']):
                        print ("Can't find ",g.demultiplexed[sample]['read1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['read2']):
                        print ("Can't find ",g.demultiplexed[sample]['read2'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index1']):
                        print ("Can't find ",g.demultiplexed[sample]['index1'])
                        exit()
                    if not os.path.isfile(g.demultiplexed[sample]['index2']):
                        print ("Can't find ",g.demultiplexed[sample]['index2'])
                        exit()

                # Bootstrap the aligned samfile paths
                # g.aligned = {}
                # for sample in g.samples:
                    # g.aligned[sample] = os.path.join(g.output_folder, 'aligned', sample + '.sam')


                g.umitag()
                g.consolidate()
                g.alignReads()
                g.identifyOfftargetSites()
                g.filterBackgroundSites()
                g.visualize()

            except Exception as e:
                print ('Error running only identify and filter.')
                print (traceback.format_exc())
                quit()
        else:
            g = GuideSeq()
            g.parseManifest(args.manifest)
            g.demultiplex()
            g.umitag()
            g.consolidate()
            g.alignReads()
            g.identifyOfftargetSites()
            g.filterBackgroundSites()
            g.visualize()

    elif args.command == 'demultiplex':
        """
        Run just the demultiplex step given the manifest
        """
        g = GuideSeq()
        g.parseManifestDemultiplex(args.manifest)
        g.demultiplex()

    elif args.command == 'umitag':
        """
        Run just the umitag step
        python guideseq/guideseq.py umitag --read1 test/data/demultiplexed/EMX1.r1.fastq --read2 test/data/demultiplexed/EMX1.r2.fastq --index1 test/data/demultiplexed/EMX1.i1.fastq --index2 test/data/demultiplexed/EMX1.i2.fastq --outfolder test/output/
        """
        g = GuideSeq()
        g.output_folder = args.outfolder
        sample = os.path.basename(args.read1).split('.')[0]
        g.samples = [sample]
        g.demultiplexed = {sample: {}}
        g.demultiplexed[sample]['read1'] = args.read1
        g.demultiplexed[sample]['read2'] = args.read2
        g.demultiplexed[sample]['index1'] = args.index1
        g.demultiplexed[sample]['index2'] = args.index2
        g.umitag()

    elif args.command == 'consolidate':
        """
        Run just the consolidate step
        python guideseq/guideseq.py consolidate --read1 test/data/umitagged/EMX1.r1.umitagged.fastq --read2 test/data/umitagged/EMX1.r2.umitagged.fastq --outfolder test/output/ --min_frequency 0.8 --min_quality 14
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.umitagged = {sample: {}}
        g.umitagged[sample]['read1'] = args.read1
        g.umitagged[sample]['read2'] = args.read2

        if 'min_quality' in args:
            min_qual = args.min_quality
        else:
            min_qual = CONSOLIDATE_MIN_QUAL

        if 'min_frequency' in args:
            min_freq = args.min_frequency
        else:
            min_freq = CONSOLIDATE_MIN_FREQ

        g.consolidate(min_freq=min_freq, min_qual=min_qual)

    elif args.command == 'align':
        """
        Run just the alignment step
        python guideseq/guideseq.py align --bwa bwa --read1 test/data/consolidated/EMX1.r1.consolidated.fastq --read2 test/data/consolidated/EMX1.r2.consolidated.fastq --genome /Volumes/Media/hg38/hg38.fa --outfolder test/output/
        """
        sample = os.path.basename(args.read1).split('.')[0]
        g = GuideSeq()
        g.BWA_path = args.bwa
        g.reference_genome = args.genome
        g.output_folder = args.outfolder
        g.samples = [sample]
        g.consolidated = {sample: {}}
        g.consolidated[sample]['read1'] = args.read1
        g.consolidated[sample]['read2'] = args.read2
        g.alignReads()

    elif args.command == 'identify':
        """
        Run just the identify step
        python guideseq/guideseq.py identify --genome /Volumes/Media/hg38/hg38.fa --aligned test/output/aligned/EMX1.sam --outfolder test/output/ --target_sequence GAGTCCGAGCAGAAGAAGAANGG
        """
        if 'description' in args:
            description = args.description
        else:
            description = ''

        if 'max_score' in args:
            max_score = args.max_score
        else:
            max_score = 7

        if 'search_radius' in args:
            search_radius = args.search_radius
        else:
            search_radius = 25

        g = GuideSeq()
        g.output_folder = args.outfolder
        g.reference_genome = args.genome
        sample = os.path.basename(args.aligned).split('.')[0]
        g.samples = {sample: {'description': description, 'target': args.target_sequence}}
        g.aligned = {sample: args.aligned}
        g.max_score = max_score
        g.search_radius = search_radius
        g.identifyOfftargetSites()

    elif args.command == 'filter':
        """
        Run just the filter step

        """
        sample = os.path.basename(args.identified).split('.')[0]
        g = GuideSeq()
        g.output_folder = args.outfolder
        g.bedtools = args.bedtools
        g.samples = {sample: {}, 'control': {}}
        g.identified = {}
        g.identified[sample] = args.identified
        g.identified['control'] = args.background
        g.filterBackgroundSites()

    elif args.command == 'visualize':
        """
        Run just the visualize step
        """
        g = GuideSeq()
        g.output_folder = os.path.dirname(args.outfolder)
        sample = os.path.basename(args.infile).split('.')[0]
        g.samples = {sample: {}}
        g.identified = {}
        g.identified[sample] = args.infile
        g.visualize()


if __name__ == '__main__':
    main()
