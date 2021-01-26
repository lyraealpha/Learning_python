#!/usr/bin/env python3.7
# -*- coding:utf-8 -*-

import os
import re
import sys
import time
import glob
import copy
import random
import requests
import argparse
import traceback
import subprocess
from bs4 import BeautifulSoup
from Ursa_Basic import timer
from Ursa_Basic import run_or_die
from Ursa_Basic import mkdir_or_die
from Ursa_Basic import capture_or_die
from Ursa_Basic import copy_to_file


def fasta_dict_parser(fasta_file, strip_id=False):
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





class KeggReceiver:

    def __init__(self, gene_id, odir='./', all_ko=None, gene_2_path=None, 
                 kegg_spename=None, pathway_anno_file=None, gene_K=None):
        self.gene_id = gene_id
        self.kegg_spename = kegg_spename
        self.odir = os.path.abspath(odir)
        self.all_ko = os.path.abspath(all_ko) if all_ko else None
        self.gene_2_path = os.path.abspath(gene_2_path) if gene_2_path else None
        self.pathway_anno_file = os.path.abspath(pathway_anno_file) if pathway_anno_file else None 
        self.gene_K = os.path.abspath(gene_K) if gene_K else None 
        mkdir_or_die(self.odir)
        self.header = {'User-Agent': 'Mozilla/5.0 (Windows NT 6.1; WOW64) AppleWebKit/537.1 (KHTML, like Gecko) '
                                'Chrome/21.0.1180.71 Safari/537.1 LBBROWSER'}

    def get_all_ko(self):
        if not self.all_ko:
            run_or_die(f'cd {self.odir} && wget -c -t 200 http://rest.kegg.jp/list/ko && '
                       f'mv ko all_ko.txt')
        if self.all_ko:
            copy_to_file(self.all_ko, f'{self.odir}/all_ko.txt')
        ko_hash = {}
        with open(f'{self.odir}/all_ko.txt', 'r') as infile:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                ko_hash[line[0]] = line[1] .lstrip(f'ko:')
        return ko_hash


    def retrieve_html_content(self, url):
        content = 'random_marker'
        while 1:
            try:
                context = requests.get(url, headers=self.header)
                content = str(BeautifulSoup(context.text, 'lxml'))
                if content !='random_marker':
                    break
            except:
                time.sleep(random.uniform(1,3))
        if content == None:
            return None
        else:
            content = content.lstrip('<html><body><p>').rstrip(f'</p></body></html>')
            return([x.split('\t')[1] for x in content.split('\n') if ':' in x])        


    def retrieve_gene_ko(self):
        url=f'http://rest.kegg.jp/link/ko/{self.kegg_spename}'
        if not self.gene_K:
            run_or_die(f'cd {self.odir} && wget -c -t 200 {url} && mv {self.kegg_spename} gene_K.txt' )
        else:
            copy_to_file(self.gene_K, f'{self.odir}/gene_K.txt')
        gene_K_hash = {}
        with open(f'{self.odir}/gene_K.txt', 'r') as infile:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                gene_K_hash[line[0]] = line[1]
        return gene_K_hash

    def get_pathways(self):

        if not self.kegg_spename:
            raise Exception('Can not get pathways if parmater kegg_spename not specified')

        if not self.gene_2_path:
            url = f'http://rest.kegg.jp/link/pathway/{self.kegg_spename}'
            run_or_die(f'cd {self.odir} && wget -c -t 200 {url} && mv {self.kegg_spename} pre_path_way.txt' )
        else:
            copy_to_file(self.gene_2_path, f'{self.odir}/pre_path_way.txt')
        pathway_hash ={}
        pth_lst = []
        with open(f'{self.odir}/pre_path_way.txt', 'r') as infile:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                gene_id = line[0]
                pth = line[1].replace(f'path:{self.kegg_spename}', 'ko')
                if not pathway_hash.get(gene_id):
                    pathway_hash[gene_id] = []
                pathway_hash[gene_id].append(pth)
                if not pth in pth_lst:
                    pth_lst.append(pth)
        if not self.pathway_anno_file:
            with open(f'{self.odir}/pre_path_anno.txt' ,'w') as ofile:
                for pth in pth_lst:
                    url = f'http://rest.kegg.jp/list/path:{pth}'
                    anno = self.retrieve_html_content(url)
                    if len(anno) != 1:
                        raise Exception(f'fail to fetch content from {url} for pathway annotation')
                    ofile.write(f'{pth}\t{anno[0]}\n')
        else:
            copy_to_file(self.pathway_anno_file, f'{self.odir}/pre_path_anno.txt')                    
        pth_anno = {}
        with open(f'{self.odir}/pre_path_anno.txt', 'r') as infile:
            for line in infile:
                line = line.rstrip('\n').split('\t')
                mapid = line[0]
                pathway_anno  = line[1]
                pth_anno[mapid] = pathway_anno
        return pathway_hash, pth_anno



def parse_kegg(pep, odir, all_ko, gene_2_path, kegg_spename, pathway_anno_file, gene_K):
    pep, odir = [os.path.abspath(x) for x in [pep, odir]]
    mkdir_or_die(odir)
    pep_d = fasta_dict_parser(pep, strip_id=1)
    gene_id = ''
    gene_pathway, pathway_anno = KeggReceiver(gene_id, odir, all_ko, gene_2_path, 
                                              kegg_spename, pathway_anno_file, gene_K).get_pathways()
    all_ko_hash = KeggReceiver(gene_id, odir, all_ko, gene_2_path, 
                               kegg_spename, pathway_anno_file, gene_K).get_all_ko()
    k_hash = KeggReceiver(gene_id, odir, all_ko, gene_2_path, 
                          kegg_spename, pathway_anno_file, gene_K).retrieve_gene_ko()
    fa_anno = f'{odir}/{kegg_spename}.fa.anno'
    fa = f'{odir}/{kegg_spename}.fa'
    pathway_ofile = f'{odir}/{kegg_spename}.gene2pathway.txt'
    n = 0
    total = len(pep_d)
    with open(f'{fa_anno}' , 'w') as ofile_fa_anno, open(fa, 'w') as ofile_fa, open(pathway_ofile, 'w') as opathway:
        for pep_id in pep_d:
            n += 1
            print(f'{round(float(n/total)*100,4)}% complete')
            if k_hash.get(pep_id):
                Knumber = k_hash[pep_id]
                if all_ko_hash.get(Knumber):
                    K_anno = all_ko_hash[Knumber]
                    Knumber = Knumber.replace(f'ko:', '')
                    anno_des = K_anno.split(';')[1].strip()
                    ofile_fa.write(f'>{pep_id} {anno_des}; {Knumber} {anno_des}\n{pep_d[pep_id]}\n')
                    ofile_fa_anno.write(f'{pep_id}\t{Knumber} {K_anno}\n')
                    if gene_pathway.get(pep_id):
                        for pathway in gene_pathway[pep_id]:
                            opathway.write(f'{pep_id}\t{pathway}\t{pathway_anno[pathway]}\n')
                    else:
                        print(f'no pathway found for gene:{pep_id}|Knumber:{Knumber}')
                else:
                    print(f'{Knumber} not found in all K file')
            else:
                print(f'{pep_id} not found in k_hash')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,
                                     description='script to get',
                                     epilog='Example: python %(prog)s')
    parser.add_argument('-pep',
                        type=str,
                        dest='pep',
                        nargs='?',
                        help='input pep file for kegg query, fasta forma; forced')
    parser.add_argument('-all_ko',
                        type=str,
                        dest='all_ko',
                        nargs='?',
                        help='all ko annotation file, automatic download if not specified; optional')
    parser.add_argument('-gene_2_path',
                        type=str,
                        dest='gene_2_path',
                        nargs='?',
                        help='gene_2_path annotation file, automatic download if not specified; optional')      
    parser.add_argument('-kegg_spename',
                        type=str,
                        dest='kegg_spename',
                        nargs='?',
                        help='kegg_spename， hsa for human; forced')        
    parser.add_argument('-pathway_anno_file',
                        type=str,
                        dest='pathway_anno_file',
                        nargs='?',
                        help='pathway_anno_file, automatic download if not specified; optional')   
    parser.add_argument('-gene_K',
                        type=str,
                        dest='gene_K',
                        nargs='?',
                        help='gene K number relation, automatic download if not specified; optional')                                                                                                                   
    parser.add_argument('-odir',
                        type=str,
                        dest='odir',
                        nargs='?',
                        default='./',
                        help='odir to for output files; default: ./')
    args = parser.parse_args()
    if not all([args.pep, args.odir]):
        parser.print_help()
        sys.exit(1)
    parse_kegg(args.pep, args.odir, args.all_ko, args.gene_2_path, args.kegg_spename, args.pathway_anno_file, args.gene_K)


