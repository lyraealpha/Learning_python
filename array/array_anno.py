#!/usr/bin/env python3.7

import os
import sys
import argparse
from Ursa_Basic import FastxParser
from Ursa_Basic import SeqManipulator
from cython_basic import smith_waterman


#SM_paramater(bc, seq, score_cutoff, identity_cutoff, bc_left_span, bc_right_span, show_match)

def create_array_anno(array_csv, genome, outfile):
    array_csv, genome, outfile = [os.path.abspath(x) for x in [array_csv, genome, outfile]]
    genome_d = FastxParser.fasta_dict_parser(genome, strip_id=1)
    with open(array_csv, 'r') as infile, open(outfile, 'w') as ofile, open(f'{outfile}.log', 'w') as _log:
        """
        assum header of csv file like below:
        IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,Species
        ,Source,SourceVersion,SourceStrand,SourceSeq,TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand
        """ 
        header =f'IlluminaID\tIlluminaName\tCHR\tPOS\tCHROM_POS\tID\tTop_Strand\tBaseOfGenome(+:/-:comp)\tSNP\tDirectionFromCsv\tsource\n'
        ofile.write(header)
        for i in infile:
            if not i.startswith('IlmnID'):
                try: 
                    (IlmnID,Name,IlmnStrand,SNP,AddressA_ID,AlleleA_ProbeSeq,
                    AddressB_ID,AlleleB_ProbeSeq,GenomeBuild,Chr,MapInfo,Ploidy,
                    Species,Source,SourceVersion,SourceStrand,SourceSeq,
                    TopGenomicSeq,BeadSetID,Exp_Clusters,RefStrand) = i.strip().split(',') 
                    if genome_d.get(Chr):
                        Ref, Alt = SNP[1:4].split('/')
                        SNX_left = TopGenomicSeq.find('[') + 25
                        SNX_right = len(TopGenomicSeq) - TopGenomicSeq.find(']') + 25
                        #extend 25 bp left and right from genome sequence make alignment more tolerant for insertion and deltion
                        csv_top_strand = TopGenomicSeq[0:TopGenomicSeq.find('[')] + TopGenomicSeq[TopGenomicSeq.find(']')+1:]
                        csv_top_strand = csv_top_strand.upper()
                        MapInfo = int(MapInfo)
                        SNX_position = MapInfo - 1
                        base_of_genome = genome_d[Chr][SNX_position]
                        if RefStrand == '-':
                            base_of_genome = SeqManipulator(base_of_genome).complement_seq()
                        SNL = SNX_position - SNX_left
                        SNR = SNX_position + SNX_right
                        ref_seq_forward = genome_d[Chr][SNL: SNR].upper()
                        ref_seq_reverse = SeqManipulator(ref_seq_forward).complement_seq(reverse=1).upper()
                        print(f'aligning reads {IlmnID}\n')
                        try:
                            print(f'align to strand +:\n')
                            al_forward = smith_waterman(csv_top_strand, ref_seq_forward, 75, 90, 30, 30, show_match=1)
                            print(f'align to strand -:\n')
                            al_reverse = smith_waterman(csv_top_strand, ref_seq_reverse, 75, 90, 30, 30, show_match=1)
                            print(f'\n\n\n')
                            Strand = 'unKnown'
                            if al_forward == 0:
                                if al_reverse == 0:
                                    _log.write(f'Can not find any match for {IlmnID}\n')
                                if al_reverse != 0:
                                    Strand = '-'
                            if al_forward != 0:
                                if al_reverse == 0:
                                    Strand = '+'
                                if al_reverse != 0:
                                    _log.write(f'Find miltiple match for {IlmnID}, unable to judge top strand\n')
                            if not Strand == 'unKnown':
                                print(f'Picking starnd : {Strand}')
                                ofile.write(f'{IlmnID}\t{Name}\t{Chr}\t{MapInfo}\t{Chr}_{MapInfo}\t.\t{Strand}\t{base_of_genome}\t{SNP}\t{RefStrand}\t{Source}\n')
                        except ZeroDivisionError:
                            print(f'Can not find any match for {IlmnID}\n')
                    else:
                        print(f'Can not find Chr_name: {Chr} in genome: {genome}') 
                except ValueError:
                    #unable to unpack if not column = len(header), then skip
                    print(f'skipping line {i}, unable to unpack')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawTextHelpFormatter,
        description='',
        epilog='Example: python %(prog)s  ')
    parser.add_argument('-array_csv',
                        type=str,
                        dest='array_csv',
                        nargs='?',
                        help='input array_csv')
    parser.add_argument('-genome',
                        type=str,
                        dest='genome',
                        nargs='?',
                        help='ref genome file') 
    parser.add_argument('-outfile',
                        type=str,
                        dest='outfile',
                        nargs='?',
                        default='./array_anno.txt',
                        help='output file, default: ./array_anno.txt')                         
    args = parser.parse_args()
    if not args.array_csv and not args.genome:
        parser.print_help()
        sys.exit(1)
    create_array_anno(args.array_csv, args.genome, args.outfile)

