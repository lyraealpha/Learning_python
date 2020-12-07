#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os
import re
import sys
import time
import argparse
import subprocess
import multiprocessing

def script_timer(defined_function):  # script timer decor
    def wrapper(*args, **kwargs):
        time_start = time.time()
        print ('\nRuning script {} starts: {}s.'.format(sys.argv[0], time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
        wrapper_out = defined_function(*args, **kwargs)
        print ('\nRuning script {} ends: {}s.'.format(sys.argv[0], time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
        time_end = time.time()
        print ('\nDone !\nRuning script {} finished: {}s elapsed.'.format(sys.argv[0], round((time_end - time_start), 4)))
        return wrapper_out
    return wrapper

def func_timer(defined_function):  # function timer decor
    def wrapper(*args, **kwargs):
        time_start = time.time()
        print ('\nRuning function {} starts: {}s.'.format(defined_function.__name__,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
        wrapper_out = defined_function(*args, **kwargs)
        print ('\nRuning function {} ends: {}s.'.format(defined_function.__name__,time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))
        time_end = time.time()
        print ('\nDone !\nRuning function {} finished: {}s elapsed.'.format(defined_function.__name__,round((time_end - time_start), 4)))
        return wrapper_out
    return wrapper

def mkdir_or_die(dir_to_make):
    target_dir = os.path.abspath(dir_to_make)
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)

def run_or_die(cmd):
    returned_ = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,shell=True)
    if returned_.returncode != 0:
        raise RuntimeError('Runing command failed ! \nFailed command:[{}]'.format(cmd))
    print('Runing command:{} finished !'.format(cmd))


def parse_orig_gff3(gff3):
    gene_dict = {}
    iso_dict = {}
    GENE_dict = {}
    gff3 = os.path.abspath(gff3)
    with open(gff3, 'r') as gff:
        f_line = gff.readlines()
        for line in f_line:
            if not line.startswith('#') and not len(line) <=5:
                line = line.strip().split('\t')
                if line[2] == 'gene':
                    gene_name = ''.join([x for x in line[8].split(';') if x.startswith('Name=')]).replace('Name=','')
                    gene_id = ''.join([x for x in line[8].split(';') if x.startswith('ID=')]).replace('ID=','')
                    gene_dict[gene_id] = '{i}({n})'.format(i=gene_id, n=gene_name)
                if line[2] == 'mRNA':
                    iso_name = ''.join([x for x in line[8].split(';') if x.startswith('Name=')]).replace('Name=','')
                    iso_id = ''.join([x for x in line[8].split(';') if x.startswith('ID=')]).replace('ID=','')
                    iso_dict[iso_id] = '{i}({n})'.format(i=iso_id, n=iso_name)

        for k,v in gene_dict.items():
            new_k = k.replace('gene','GENE')
            new_v = v.replace('gene','GENE')
            GENE_dict.update({k: v})
            GENE_dict.update({new_k:new_v})
    return GENE_dict,iso_dict


def reform_regrex_pattern_match_list(in_list):
    final_list = []
    if in_list[0].startswith('rna'):
        tmp_list = sorted([int(x.replace('rna','')) for x in in_list],reverse=True)
        for i in [''.join(['rna',str(x)]) for x in tmp_list]:
            if not i in final_list:
                final_list.append(i)
    if in_list[0].startswith('gene'):
        tmp_list = sorted([int(x.replace('gene','')) for x in in_list],reverse=True)
        for i in [''.join(['gene',str(x)]) for x in tmp_list]:
            if not i in final_list:
                final_list.append(i)
    if in_list[0].startswith('GENE'):
        tmp_list = sorted([int(x.replace('GENE','')) for x in in_list],reverse=True)
        for i in [''.join(['GENE',str(x)]) for x in tmp_list]:
            if not i in final_list:
                final_list.append(i)
    return final_list

def run_pattern_match_do_relpacement_else_return(regrex_pattern,in_line,in_dict):
    template_line = in_line
    line_replacable = []
    find_target = regrex_pattern.findall(in_line)
    if len(find_target) > 0:
        reordered_find_target = reform_regrex_pattern_match_list(find_target)
        for target in reordered_find_target:
            if in_dict.get(target):
                line_replacable.append('True')
                template_line = template_line.replace(target,in_dict[target])
    replace_bool = True if 'True' in line_replacable else False
    return template_line,replace_bool

def file_line_processing(in_line,gene_dict,rna_dict):
    GENE_regrex = re.compile(r'GENE\d+')
    gene_regrex = re.compile(r'gene\d+')
    rna_regrex = re.compile(r'rna\d+')
    (GENE_line, GENE_replacable_bool) = run_pattern_match_do_relpacement_else_return(GENE_regrex,in_line,gene_dict) # rename gene first
    (gene_line, gene_replacable_bool) = run_pattern_match_do_relpacement_else_return(gene_regrex, GENE_line[:], gene_dict)  # rename gene first
    (rna_line, rna_replacable_bool) = run_pattern_match_do_relpacement_else_return(rna_regrex, gene_line[:], rna_dict) # then rename rna followed
    final_line = rna_line[:]
    final_bool = True if any(([GENE_replacable_bool ==True],[gene_replacable_bool ==True],[rna_replacable_bool ==True])) else False
    return final_line,final_bool

def file_processing(infile,gene_dict,rna_dict,id,od):
    infile = os.path.abspath(infile)
    oufile = infile.replace(id, od)
    (oufile_dir, oufile_name) = os.path.split(oufile)
    if not os.path.isdir(oufile_dir):
        os.makedirs(oufile_dir)
    only_copy_condition = [(infile.endswith('png')),
                           (infile.endswith('pdf')),
                           (infile.endswith('.gz')),
                           (infile.endswith('.css')),
                           (infile.endswith('.less')),
                           ('js/plugins' in infile),
                           ('Heatmap_HTML' in infile),
                           ('/js' in infile),
                           (infile.endswith('.scss')),
                           (infile.endswith('.css.map')),
                           (infile.endswith('.otf')),
                           (infile.endswith('.xml')),
                           (infile.endswith('.eot')),
                           (infile.endswith('.ttf')),
                           (infile.endswith('.woff')),
                           (infile.endswith('.woff2')),
                           (infile.endswith('HTML/All.DEG_final')),
                           (infile.endswith('.js')),
                           (infile.endswith('.svg')),
                           (infile.endswith('All_Database_annotation.xls')),
                           (infile.endswith('.jpg')),
                           (infile.endswith('.gif')),
                           (infile.endswith('.cloud')),
                           (infile.endswith('html') and infile.startswith('ko'))]
    if any(only_copy_condition):
        run_or_die('cp {} {}'.format(infile,oufile))
        print('Only copy file:{} found,just copy !'.format(infile))
        return
    replaceable = []
    fixed_file = []
    with open(infile, 'r') as input:
        f_line = input.readlines()
        for line in f_line:
            (fixed_line,replace_bool) = file_line_processing(line,gene_dict,rna_dict)
            replaceable.append(replace_bool) if replace_bool == 'True' else None
            fixed_file.append(fixed_line)
    with open(oufile, 'w') as out:
        out.write(''.join(fixed_file))

def walk_through_dir_and_return_all_file(id):
    all_file_list = []
    in_dir = os.path.abspath(id)
    three_dimen_tuple_list = os.walk(in_dir)
    for item in three_dimen_tuple_list:
        (root_dir,sub_dirs,files) = item
        if len(files) > 0:
            [all_file_list.append(os.path.join(root_dir,file)) for file in files]
    return all_file_list


@script_timer
def generate_result(indir,outdir,gff3):
    in_dir, out_dir, gff3 = os.path.abspath(indir), os.path.abspath(outdir), os.path.abspath(gff3)
    (gene_dict,rna_dict) = parse_orig_gff3(gff3)
    all_file_list = walk_through_dir_and_return_all_file(in_dir)
    tmp_list = []
    pool = multiprocessing.Pool(10)
    for file in all_file_list:
        tmp_list.append(pool.apply_async(file_processing,args=(file,gene_dict,rna_dict,in_dir,out_dir)))
    pool.close()
    pool.join()





if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='Walk through dir, and add gene name column to files if needed ',
        epilog='Example: python %(prog)s -id  -od   -gff3  species_name.gff3')
    parser.add_argument('-id',
                        type=str,
                        dest='indir',
                        nargs='?',
                        help='input dir (source)')
    parser.add_argument('-od',
                        type=str,
                        dest='outdir',
                        nargs='?',
                        help='modified result dir')
    parser.add_argument('-gff3',
                        type=str,
                        dest='gff3',
                        nargs='?',
                        help='gff3 file , gene_id and gene_name are both required in feature line')
    args = parser.parse_args()
    if not all([args.indir,args.outdir,args.gff3]):
        parser.print_help()
        sys.exit(1)
    generate_result(args.indir, args.outdir, args.gff3)
