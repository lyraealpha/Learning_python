#!/usr/bin/env python3.7
# -*- coding:utf-8 -*-

import os
import re
import sys
import time
import glob
import copy
import requests
import argparse
import traceback
import subprocess
from bs4 import BeautifulSoup
from Ursa_Basic import timer
from Ursa_Basic import run_or_die
from Ursa_Basic import mkdir_or_die
from Ursa_Basic import capture_or_die


def download_or_die(cmd):  # python 3.6+
    time_start = time.time()
    run_cmd = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    return_code = run_cmd.returncode
    if return_code != 0:
        errors = run_cmd.stderr.decode('utf-8')
        raise DownLoadingSraFileError(f'Running command failed ! \n'
                                      f'Failed command: {cmd}\n'
                                      f'Error description:\n'
                                      f'{errors}')
    time_end = time.time()
    print('Running command : {} finished: {}s elapsed.'.format(cmd, round((time_end - time_start), 4)), flush=True)


class DownLoadingSraFileError(Exception):

    def __init__(self, e='Error while Downloading SRA file!'):
        self.e = e
        super().__init__()

    def __str__(self):
        return self.e


class SraReceiver:

    def __init__(self, sra_name, argument_l, tmp_dir):
        self.tmp_dir = os.path.abspath(tmp_dir)
        self.argument_l = argument_l
        self.sra_name = sra_name
        self.ascp = '/opt/aspera/connect/bin/ascp'
        self.ascp_key = '/opt/aspera/connect/etc/asperaweb_id_dsa.openssh'
        self.target_ftp = (f'anonftp@ftp-trace.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/'
                           f'sra/{sra_name[0:3]}/{sra_name[0:6]}/{sra_name}/{sra_name}.sra')

    def fasp_method_receiver(self):
        mkdir_or_die(self.tmp_dir，clean_before_make=1)
        ascp_cmd = (f'{self.ascp} '
                    f'-T '
                    f'-k 1 '
                    f'-i {self.ascp_key} '
                    f'-l {self.argument_l} {self.target_ftp} {self.tmp_dir} ')
        print(f'try to download SRA file (ascp) for {self.sra_name}')
        for try_time in [x for x in range(1, 11)]:
            try:
                # print(f'Try to download SRA file {self.sra_name} via ascp in round{try_time}')
                download_or_die(ascp_cmd)
                return 'Succeed'
            except DownLoadingSraFileError:
                pass
        return 'Failed'

    @staticmethod
    def target_sra_https_finder(sra_name):  # get SRA download link from ncbi
        target_url = f'https://trace.ncbi.nlm.nih.gov/Traces/sra/?run={sra_name}'
        header = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) '
                                'Chrome/21.0.1180.71 Safari/537.1 LBBROWSER'}
        context = requests.get(target_url, headers=header)
        content = BeautifulSoup(context.text, 'lxml')
        download_url = content.find_all(href=re.compile(f'sra-download.*{sra_name}'))
        resources = str(download_url).split('/a>')
        sra_links = []
        for resource in resources:
            sra_link = re.match('.*"(.*)".*>?', str(resource))
            if sra_link:
                sra_links.append(sra_link.group(1))
        if len(sra_links) < 1:
            raise Exception(f'Can not find SRA file download link for: {sra_name} ')
        return sra_links

    def wget_method_receiver(self):  # try wget way
        clean_and_mkdir(self.tmp_dir)
        sra_links = SraReceiver.target_sra_https_finder(self.sra_name)
        for sra_link in sra_links:
            print(f'try to download SRA file (wget) for {self.sra_name} from {sra_link}')
            round_start = f'cd {self.tmp_dir} && wget --tries=200 {sra_link}'
            try:
                download_or_die(round_start)
                return 'Succeed'
            except DownLoadingSraFileError:
                pass
            for try_time in [x for x in range(1, 6)]:
                # print(f'Try to download SRA file {self.sra_name} via wget in round{try_time}')
                round_end = f'cd {self.tmp_dir} && wget --continue --tries=200 {sra_link}'
                try:
                    download_or_die(round_end)
                    return 'Succeed'
                except DownLoadingSraFileError:
                    pass
        return 'Failed'

    @timer
    def sra_receiver(self):
        try_ascp = SraReceiver(self.sra_name, self.argument_l, self.tmp_dir).fasp_method_receiver()
        if try_ascp == 'Succeed':
            return 'Succeed'
        elif try_ascp == 'Failed':
            try_wget = SraReceiver(self.sra_name, self.argument_l, self.tmp_dir).wget_method_receiver()
            if try_wget == 'Succeed':
                return 'Succeed'
            if try_wget == 'Failed':
                raise DownLoadingSraFileError


def receiver(sra_name, argument_l, tmp_dir):
    tmp_dir = os.path.abspath(f'{tmp_dir}/{sra_name}')
    SraReceiver(sra_name, argument_l, tmp_dir).sra_receiver()
    download_file = glob.glob(f'{tmp_dir}/*{sra_name}*')[0]
    if not download_file.endswith('.sra'):
        run_or_die(f'mv {download_file} {tmp_dir}/{sra_name}.sra')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='script to down load SRA archive file from NCBI',
                                     epilog='Example: python %(prog)s -sra_name SRR1649426 '
                                            '-argument_l 5M -tmp_dir /path/to/target/dir')
    parser.add_argument('-sra_name',
                        type=str,
                        dest='sra_name',
                        nargs='?',
                        help='input sra_name; SRR1649426 or ERR2984736')
    parser.add_argument('-argument_l',
                        type=str,
                        dest='argument_l',
                        nargs='?',
                        default='10M',
                        help='argument_l for ascp download method')
    parser.add_argument('-tmp_dir',
                        type=str,
                        dest='tmp_dir',
                        nargs='?',
                        help='tmp_dir for download SRA file , download file in tmpdir/sra_name')
    args = parser.parse_args()
    if not all([args.sra_name, args.argument_l, args.tmp_dir]):
        parser.print_help()
        sys.exit(1)
    receiver(args.sra_name, args.argument_l, args.tmp_dir)
