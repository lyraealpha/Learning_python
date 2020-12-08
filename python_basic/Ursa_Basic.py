#!/usr/bin/env python3.7
# -*- coding:utf-8 -*-

"""
Basic Functions 
"""

import os
import re
import sys
import time
import glob
import copy
import math
import random
import argparse
import traceback
import subprocess
import multiprocessing
from cython_basic import qv_calculator

try:  # pb
    import xlrd
    import pysam
except ModuleNotFoundError:
    pass


def find_hifidelity(ifq, ofq):
    ifq = os.path.abspath(ifq)
    ofq = os.path.abspath(ofq)
    with open(ifq, 'r') as infile, open(ofq, 'w') as ofile:
        while 1:
            seq_id = infile.readline()
            if seq_id == '':
                break
            seq, sup_info, qvs = infile.readline(), infile.readline(), infile.readline()
            predict_accuracy, qv = qv_calculator(qvs.strip())
            if predict_accuracy > 0.99:
                ofile.write(f'{seq_id.rstrip()}\tacc={predict_accuracy}\tqv={qv}\n{seq}{sup_info}{qvs}')


def timer(defined_function):  # time decorator
    def wrapper(*args, **kwargs):
        time_start = time.time()
        wrapper_out = defined_function(*args, **kwargs)
        time_end = time.time()
        print('Done !\nRuning script {} finished: {}s elapsed.'.format(sys.argv[0], round((time_end - time_start), 4)))
        return wrapper_out

    return wrapper


def find_abs_path(*args):
    lst = []
    for i in args:
        try:
            i = os.path.abspath(i) if os.path.isfile(i) or os.path.isdir(i) else i
        except:
            pass
        lst.append(i)
    return (lst)


class RunOrDieError(Exception):
    def __init__(self, e='Running command failed'):
        self.e = e
        super().__init__()

    def __str__(self):
        return self.e


def run_or_die(cmd, time_track=True):
    time_start = time.time()
    run_cmd = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return_code = run_cmd.returncode
    if return_code != 0:
        errors = run_cmd.stderr.decode('utf-8')
        raise RunOrDieError(f'Running command failed ! \n'
                            f'Failed command: {cmd}\n'
                            f'Error description:\n'
                            f'{errors}')
    if time_track:
        time_end = time.time()
        print('Running command : {} finished: {}s elapsed.'.format(cmd, round((time_end - time_start), 4)), flush=True)

def mkdir_or_die(dir_to_make, clean_before_make=False):
    target_dir = os.path.abspath(dir_to_make)
    if len(target_dir) < 15:
        raise RuntimeError(f'Target dir too short(at least longer than '
                            f'15,'
                            f'stop in case unexcepted error')
    if os.path.isdir(target_dir) and clean_before_make:
        run_or_die('rm -r {}'.format(target_dir))
    try:
        os.makedirs(target_dir)
    except FileExistsError:  # in case of Race Condition
        pass

def mkdirs_or_dies(*args):
    for i in args:
        try:
            mkdir_or_die(str(i))
        except Exception as e:
            print(f'creating directory failed {i}: \nError:{e}')

class CaptureFailedError(Exception):

    def __init__(self, e='Progress failed in capture results from command'):
        self.e = e
        super().__init__()

    def __str__(self):
        return self.e

def capture_or_die(cmd):
    """
    capture result form shell command
    Raise exception and return failed command if any error occurred
    Note: return raw results, '\n' not striped or processed
    """
    run_cmd = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return_result = run_cmd.stdout
    return_code = run_cmd.returncode
    if return_code == 0:
        return return_result.decode('utf-8')
    if return_code != 0:
        errors = run_cmd.stderr.decode('utf-8')
        raise CaptureFailedError(f'Running command failed! \n'
                                    f'Failed command: {cmd}!\n'
                                    f'Error description: \n'
                                    f'{errors}\n')

def capture_errors(cmd, verbose=False, raise_if_failed=True):
    """
    capture errors form shell command
    Raise exception and return failed command if any error occurred
    mainly used in multiprocessing child process error handling
    """
    run_cmd = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return_code = run_cmd.returncode
    if raise_if_failed == True:
        if return_code != 0:
            errors = run_cmd.stderr.decode('utf-8')
            if not verbose:
                raise CaptureFailedError(f'Running command failed! \n Failed command: {cmd}!\n'
                                            f'Error description: \n {errors}\n')
            else:
                raise CaptureFailedError(f'Running command failed! \n Failed command: {cmd}!\n'
                                            f'Error description: \n {errors}\n'
                                            f'{traceback.format_exc()}\n')

def copy_to_dir(resource, target):
    """
    resource can be a list of files/directories; or can be one file/direcrory
    """
    target = os.path.abspath(target)
    source = glob.glob(f'{resource}') if not isinstance(resource, list) else resource
    if not os.path.isdir(target):
        mkdir_or_die(target)
    if len(source) == 1:
        for item in source:
            item = os.path.abspath(item)
            # incase parsing wrong directory (target dir shorter than 15 will raise Exception)
            if len(item) < 15:
                raise Exception(f'Target dir too short(at least longer than '
                                f'15,'
                                f'stop in case unexcepted error')
            if os.path.isdir(item):
                run_or_die(f'cp -r {item} {target}')
            if os.path.isfile(item):
                run_or_die(f'cp {item} {target}')
    if len(source) > 1:
        for item in source:
            item = os.path.abspath(item)
            if len(item) < 15:
                raise Exception(f'Target dir too short(at least longer than '
                                f'15,'
                                f'stop in case unexcepted error')
            if os.path.isdir(item):
                run_or_die(f'cp -r {item} {target}')
            if os.path.isfile(item):
                run_or_die(f'cp {item} {target}')
    if len(source) == 0:
        print(f'No resource found for :{resource}, skip copy')

def copy_to_file(resource, target):
    target = os.path.abspath(target)
    source = glob.glob(f'{resource}')
    if len(source) > 1:
        raise Exception('Unable to copy multi files to one file')
    if len(source) == 1:
        target_base_dir = os.path.split(target)[0]  # incase target file dir not exists
        if not os.path.isdir(target_base_dir):
            mkdir_or_die(target_base_dir)
        run_or_die(f'cp {resource} {target}')
    if len(source) == 0:
        print(f'No resource found for :{resource}, skip copy')

def runs_or_dies(object_, threads=10, verbose=False):
    """
    :param object_: list[run_cmd1,run_cmd2...] or shell file(run_cmds.sh)
    :param threads: threads for multiprocessing
    :param verbose: error details
    :return: None
    """
    time_start = time.time()
    pool = multiprocessing.Pool(threads)
    runs = []
    if isinstance(object_, list):
        command_list = copy.deepcopy(object_)
        for cmd in command_list:
            runs.append(pool.apply_async(capture_errors, args=(cmd, verbose)))
        pool.close()
        pool.join()
    else:
        if os.path.isfile(object_):
            with open(os.path.abspath(object_), 'r') as infile:
                for cmd in [line.strip() for line in infile.readlines() if len(line) >= 2]:
                    runs.append(pool.apply_async(capture_errors, args=(cmd, verbose)))
            pool.close()
            pool.join()
    checker = f'marker:{time_start}'
    for child_res in runs:
        checker = child_res
        try:
            checker.get()
        except Exception as e:
            if not verbose:
                raise Exception(f'{e}\n')
            else:
                raise Exception(f'{e}\n{traceback.format_exc()}')
    if checker == f'marker:{time_start}':
        raise Exception(f'No child res, wrong object_ type or commands:\n'
                        f'{object_}')
    print(f'Commands done, time elapsed: {round((time.time() - time_start), 4)} s')

class QsubHandler:

    def __init__(self):
        pass

    @staticmethod
    def qsub_or_die(shell_cmd, queue, cpu_num, resource_vf, check_qsub=False):
        if not shell_cmd:
            raise Exception('Undefined shell file {}'.format(shell_cmd))
        sh_cmd = os.path.abspath(shell_cmd)
        if not 1 <= int(cpu_num) <= 100:
            raise Exception('Wrong cpu number defined {}, -maxproc should in [1 : 100]'.format(cpu_num))
        init_cmd = 'sh path/to/qsub_sge.plus.sh  '
        init_cmd += '--queue {q} --maxproc {c} --independent --reqsub  '.format(q=queue, c=cpu_num)
        if resource_vf:
            init_cmd += '--resource vf={v}  '.format(v=resource_vf)
        init_cmd += f'{sh_cmd}'
        print(init_cmd)
        run_or_die(init_cmd)
        if check_qsub:
            QsubHandler.qsub_check(shell_cmd)

    @staticmethod
    def qsub_check(shfile):
        sh_file = os.path.abspath(shfile)
        qusub_sh = glob.glob('{}*.qsub/*.sh'.format(sh_file))
        checked_sh = glob.glob('{}*.qsub/*.Check'.format(sh_file))
        for _sh in qusub_sh:
            if not ''.join([_sh, '.Check']) in checked_sh:
                raise Exception('Warning ! qsub check failed in: {}'.format(_sh))

class FastxParser:

    def __init__(self, fastx_file):
        self.fastx_file = os.path.abspath(fastx_file)
        self.f_type = 'fastq' if any([fastx_file.endswith('fastq'), fastx_file.endswith('fq')]) else 'fasta'

    @staticmethod
    def fasta_dict_parser(fasta_file, strip_id=False):
        """
        :param fasta_file: input fasta format file
        :param strip_id: strip redundant information from  (>chr1 xxx-xxx xxx xxxx) to (chr1)
        :return: {id : seq} dictionary
        |incase not type below                           |incase type below
        |>seq1                                           |>seq1
        |AGCTN                                           |AGCTN
        |>seq2                                           |NAGCT
        |NAGCT...                                        |>seq2...
        """
        fa_dict = {}
        seq_id = ''
        with open(fasta_file, 'r') as infastx:
            while 1:
                f_line = infastx.readline()
                if f_line == '':
                    break
                if f_line.startswith('>'):
                    if strip_id:
                        seq_id = f_line.strip().split()[0][1:]
                    else:
                        seq_id = f_line
                    fa_dict[seq_id] = []
                else:
                    fa_dict[seq_id].extend([f_line.strip()])
        for seq_id, sequence in fa_dict.items():
            fa_dict[seq_id] = ''.join(sequence)
        return fa_dict

    def base_summer(self):
        base_num = 0
        seq_num = 0
        if self.f_type == 'fastq':
            with open(self.fastx_file, 'r') as infastx:
                while 1:
                    read_id = infastx.readline()
                    if read_id == '':
                        break
                    if not read_id.startswith('@'):
                        raise Exception(f'Wrong line file for fastq read id {read_id}')
                    sequence = infastx.readline().strip()
                    base_num += len(sequence)
                    seq_num += 1
                    infastx.readline()  # +
                    infastx.readline()  # quality score line
        elif self.f_type == 'fasta':
            with open(self.fastx_file, 'r') as infastx:
                while 1:
                    f_line = infastx.readline()
                    if f_line == '':
                        break
                    if f_line.startswith('>'):
                        seq_num += 1
                    elif not f_line.startswith('>'):
                        base_num += len(f_line.strip())
        average_seq_len = base_num / seq_num
        return base_num, seq_num, average_seq_len


class SeqManipulator:  # sequences handler
    def __init__(self, in_sequence):
        self.in_sequence = in_sequence

    def complement_seq(self, reverse=False):
        f_ = 'agctAGCTnN'
        r_ = 'tcgaTCGAnN'
        dict_ = dict(zip(f_, r_))
        if reverse:
            return ''.join(map(lambda x: dict_[x], self.in_sequence))[::-1]
        else:
            return ''.join(map(lambda x: dict_[x], self.in_sequence))

    def base_percent_content(self):
        seq_len = len(self.in_sequence)
        fortmat_str = ''
        for i in 'AGCTN':
            count_x = seq.count(base.upper()) + seq.count(base.lower())
            fortmat_str += (f'{base.upper()}:{count_x} ; 
                                percent: {round(float(count_x / len(seq)) * 100, 2)}% \n')
        return fortmat_str

    @staticmethod
    def calculate_predict_accuracy(in_sequence):
        error_sum = 0
        for char in in_sequence:
            q = ord(char)
            error_possibility = 10 ** (float(q - 33) / -10)
            error_sum += error_possibility
        predict_accuracy = 1 - error_sum / len(in_sequence)
        return predict_accuracy

    @staticmethod
    def calculate_average_qscore(in_string):
        return -10 * math.log(sum([10 ** (float((ord(j)) - 33) / -10) for j in in_string]) / len(in_string), 10)

class XamParser:

    def __init__(self, in_bam_or_sam_file):
        self.infile = os.path.abspath(in_bam_or_sam_file)
        self.in_bam_file = pysam.AlignmentFile(self.infile, 'rb', check_sq=False)

    @staticmethod
    def xam_record_profiler(show_code2_str=False):
        """
        0 -- (possioble of being mapped to ref strand)
        flag description from :
        https://github.com/Magdoll/cDNA_Cupcake/blob/c2bb4b8425f5950792d75083ee3a4eed7e75f527/sequence/BioReaders.py
        Heng Li's SAM https://samtools.github.io/hts-specs/SAMv1.pdf
        1 -- 0x1 -- read is one of a pair
        2 -- 0x2 -- alignment is one end of proper PE alignment          (IGNORE)
        4 -- 0x4 -- read has no reported alignments                      (IGNORE)
        8 -- 0x8 -- read is one of a pair and has no reported alignments (IGNORE)
        16 -- 0x10 -- reverse ref strand
        32 -- 0x20 -- other mate is aligned to ref strand
        64 -- 0x40 -- first mate in pair
        128 -- 0x80 -- second mate in pair
        256 -- 0x100 -- not primary alignment
        512 -- 0x200 -- not passing QC filters
        1024 -- 0x400 -- PCR or optical duplicate
        2048 -- 0x800 -- supplementary alignment

        F0X900 -- remove supplementary or secondary alignment records
        code2str from https://github.com/rasto2211/nanopore-read-align/blob/master/scripts/sam_utils.py
        """
        code2str = {'M': "MATCH",  # M - alignment match (can be a sequence match or mismatch)
                    'I': "INS",  # I - insertion to the reference
                    'D': "DEL",  # D - deletion from the reference
                    'N': "SKIP",  # N - skipped region from the reference
                    'S': "SOFT_CLIP",  # S - soft clipping (clipped sequences present in SEQ)
                    'H': "HARD_CLIP",  # H - hard clipping (clipped sequences NOT present in SEQ)
                    'P': "PAD",  # P - padding (silent deletion from padded reference)
                    '=': "EQUAL",  # = - sequence match
                    'X': "DIFF",  # X - sequence mismatch
                    }
        if show_code2_str:
            for k, v in code2str.items():
                print(f'{k}\t{v}')

    def extract_zmws(self, extract_list, odir):
        """
        extract_list format
        #Sample_id\tZMW_id(movie included)\n
        Sample_name\tm64054_200710_131319/0\n
        """
        odir = os.path.abspath(odir)
        mkdir_or_die(odir)
        sample_zmws_dict = {}
        file_d = {}
        with open(extract_list, 'r') as sample_2_zmws:
            for line in sample_2_zmws:
                if not line.startswith('#'):
                    sample, zmw_id = line.strip().split()
                    if not sample_zmws_dict.get(zmw_id):
                        sample_zmws_dict[zmw_id] = sample
                    if not file_d.get(sample):
                        file_d[sample] = pysam.AlignmentFile(f'{odir}/{sample}.subreads.bam', 'wb',
                                                                template=self.in_bam_file)
        for reads in self.in_bam_file:
            zmw_name = '/'.join(reads.qname.split('/')[:2])
            if sample_zmws_dict.get(zmw_name):
                file_d[sample_zmws_dict[zmw_name]].write(reads)
        [file_d[f].close() for f in file_d]

    def subreads_number_2_zmw_dict(
            self):  # return zmw_number --> subreads numbers; SMRT_LINK_7.0 reclaim singletons
        find_num_passes = re.compile(r'\(\'np\'.*?\)')
        zmw_subreads_number_dict = {}
        for records in self.in_bam_file:
            zmw_id = str(records).split()[0]
            num_full_passes = find_num_passes.findall(str(records))
            if 'np' not in num_full_passes[0]:
                raise Exception('wrong np found {};failed to find full pass number for ID:{} in file :{}'
                                .format(num_full_passes, zmw_id, self.infile))
            if not len(num_full_passes) == 1:
                raise Exception('wrong number of np found {};failed to find full pass number for ID:{} in file :{}'
                                .format(num_full_passes, zmw_id, self.infile))
            np = int(''.join(num_full_passes[0].lstrip('\"(\'np\', ').rstrip(')')))
            if np >= 6:
                zmw_subreads_number_dict[zmw_id] = np
        self.in_bam_file.close()
        return zmw_subreads_number_dict

def dir_walker(idir, file_prefix=None, file_suffix=None):
    idir = os.path.abspath(idir)
    file_list = []
    dir_items = os.walk(idir)
    for dir_item in dir_items:
        root_dir, _, files = dir_item
        if len(files) > 0:
            if not file_prefix and not file_suffix:
                [file_list.append(f'{root_dir}/{f_name}') for f_name in files]
            if file_prefix:
                [file_list.append(f'{root_dir}/{f_name}') for f_name in files if f_name.startswith(file_prefix)]
            if file_suffix:
                [file_list.append(f'{root_dir}/{f_name}') for f_name in files
                    if f_name.endswith(file_suffix) and f'{root_dir}/{f_name}' not in file_suffix]
    return file_list

class SpeedyLineCollector:
    """
    infile: input file (text file only)
    threads: threads for multiprocessing
    pointer_unit: chunk size for jump reading
    collect_pointer: rerurn pointers if True
    Note: compressed file not supported
    """

    def __init__(self, infile, threads=10, collect_pointer=None):
        self.infile = os.path.abspath(infile)
        self.threads = threads if isinstance(threads, int) else int(threads)
        self.collect_pointer = collect_pointer if collect_pointer else None

    @staticmethod
    def file_reader(in_file, pointer, pointer_unit, eof):
        pointer_lst = []
        chunk_end = pointer + pointer_unit
        with open(in_file, 'r') as infile:
            infile.seek(pointer)
            if pointer != 0:
                """
                other chunks (except 1st chunk)
                skip end line of chunck before
                try to keep last pointer complete, read one more line. 
                and reading for chunk will skip the first line(complete or incomplete) to adjust cond1 and cond2 both.
                """
                infile.readline()
            while 1:
                cur = infile.tell()
                infile.readline()
                if cur > chunk_end or cur == eof:
                    """
                    cond 1: if cur = chunk_end: means the ending_point(pointer+pointer_unit) locate exactly at end of
                            current chunk and at the begining of the next chunk, 
                            pointer for line reading in next chunk will be complete, no need to skip
                    cond 2: if cur > chunk_end: means (pointer+pointer_unit) inside last line for a chunk, 
                            pointer for line reading in next chunk will be incomplete , need to skip
                    """
                    print('End of chunk reached, break')
                    break
                pointer_lst.append(cur)
        return pointer_lst

    @staticmethod
    def estimate_eof(infile):
        infile = os.path.abspath(infile)
        f_size = os.path.getsize(infile)
        return f_size

    def collects_or_dies(self):
        pool = multiprocessing.Pool(self.threads)
        runs = []
        end_point = SpeedyLineCollector.estimate_eof(self.infile)
        pointer_unit = int(end_point / 10) + 1
        chunk_lst = [x for x in range(0, end_point, pointer_unit)]
        collected_pointers = []
        time_start = time.time()
        for chunk in chunk_lst:
            runs.append(
                pool.apply_async(SpeedyLineCollector.file_reader,
                                    args=(self.infile, chunk, pointer_unit, end_point)))
        pool.close()
        pool.join()
        checker = f'marker:{time_start}'
        for child_res in runs:
            checker = child_res
            try:
                collected_pointers.extend(checker.get())
            except Exception as e:
                raise Exception(f'{e}\n')
        if checker == f'marker:{time_start}':
            raise Exception(f'No child res while processing file {self.infile}\n')
        print(
            f'Collecting pointers from file {self.infile} done, time elapsed: {round((time.time() - time_start), 4)} s\n')
        return collected_pointers

class GffParser:
    def __init__(self, gff_file):
        (self.file_path, self.file_name) = os.path.split(os.path.abspath(gff_file))
        self.gff_file = os.path.abspath(gff_file)

    def gff_reader(self):
        ID_regx = re.compile(r'ID=(.*)')
        Parent_regx = re.compile(r'Parent=(.*)')

        def find_feature_id(features, regrex):
            found_id = []
            for i in features.split(';'):
                j = regrex.match(i) if regrex.match(i) else None
                if j:
                    found_id.append(j.group(1))
            if len(found_id) == 0:
                raise Exception(f'Can not find match regx id in item {features}')
            if len(found_id) > 1:
                raise Exception(f'Find more than one match regx id in item {features}')
            return found_id[0]

        def get_pos(lst):
            """
            chr strand  start   end
            """
            return [lst[0], lst[6], int(lst[3]) - 1, int(lst[4])]

        gene_pos, mRNA_pos, cds_pos, exon_pos = {}, {}, {}, {}

        with open(self.gff_file, 'r') as infile:
            for line in infile:  # [0:chr  1:source 2:type 3:start 4:end 5:. 6:strand 7:? 8:feature]
                if not line.startswith('#'):
                    line = line.strip().split('\t')
                    if len(line) == 9:
                        if int(line[3]) > int(line[4]):
                            raise Exception('corruped line in gff_file : {}, start > end'.format(line))
                        if line[2] == 'gene':  # gene_id ----> [chr, strand, Start, End]
                            gene_id = find_feature_id(line[8], ID_regx)
                            gene_pos[gene_id] = get_pos(line)
                        if line[2] == 'mRNA':  # mRNA_id ----> [chr, strand, Start, End]
                            mRNA_id = find_feature_id(line[8], ID_regx)
                            mRNA_pos[mRNA_id] = get_pos(line)
                        if line[2] == 'exon':  # mRNA_id ----> [[chr, strand, Start, End]......]
                            exon_parent_id = find_feature_id(line[8], Parent_regx)
                            if not exon_pos.get(exon_parent_id):
                                exon_pos[exon_parent_id] = []
                            exon_pos[exon_parent_id].append(get_pos(line))
                        if line[2] == 'CDS':  # mRNA_id ----> [[chr, strand, Start, End]......]
                            cds_parent_id = find_feature_id(line[8], Parent_regx)
                            if not cds_pos.get(cds_parent_id):
                                cds_pos[cds_parent_id] = []
                            cds_pos[cds_parent_id].append(get_pos(line))
        return gene_pos, mRNA_pos, cds_pos, exon_pos

    def get_seq_dicts(self, genome_file):
        (gene_pos, mRNA_pos, cds_pos, exon_pos) = self.gff_reader()
        f_d = FastxParser.fasta_dict_parser(genome_file, strip_id=1)
        gene, mRNA, cds, exon = {}, {}, {}, {}

        def fetch_sequence(fasta_dict, lst):
            chromosome, _, start, end = lst
            sequence = fasta_dict[chromosome]
            return sequence[start: end]

        def reform_array(lst):
            return sorted(lst, key=lambda x: x[2])

        for k, v in gene_pos.items():  # gene
            if v[1] == '+':
                gene[k] = fetch_sequence(f_d, v)
            else:
                gene[k] = SeqManipulator(fetch_sequence(f_d, v)).complement_seq(reverse=1)
        for k, v in mRNA_pos.items():  # mRNA
            if v[1] == '+':
                mRNA[k] = fetch_sequence(f_d, v)
            else:
                mRNA[k] = SeqManipulator(fetch_sequence(f_d, v)).complement_seq(reverse=1)
        for k, v in cds_pos.items():  # cds
            if v[0][1] == '+':
                cds[k] = ''.join([fetch_sequence(f_d, _cds) for _cds in v])
            else:
                v = reform_array(v)
                cds[k] = SeqManipulator(''.join([fetch_sequence(f_d, _cds) for _cds in v])).complement_seq(
                    reverse=1)
        for k, v in exon_pos.items():  # exon
            if v[0][1] == '+':
                exon[k] = ''.join([fetch_sequence(f_d, _exon) for _exon in v])
            else:
                v = reform_array(v)
                exon[k] = SeqManipulator(''.join([fetch_sequence(f_d, _exon) for _exon in v])).complement_seq(
                    reverse=1)
        return gene, mRNA, cds, exon

    def extract_seqeuence(self, genome_file, odir):

        def write_fasta(ofile_name, seq_dict):
            with open(ofile_name, 'w') as ofile:
                for k, v in seq_dict.items():
                    ofile.write(f'>{k}\n{v}\n')

        mkdir_or_die(odir)
        gene, mRNA, cds, exon = self.get_seq_dicts(genome_file)
        write_fasta(f'{odir}/gene.fasta', gene)
        write_fasta(f'{odir}/mRNA.fasta', mRNA)
        write_fasta(f'{odir}/cds.fasta', cds)
        write_fasta(f'{odir}/exon.fasta', exon)

    if __name__ == "__main__":
        pass
