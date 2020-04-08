#! /usr/bin/env python3
#This script searches for an open reading frame and translates it.
#Call on script: python translate_orf.py CGUACGUAACGUACUGCAGUCAUACGUACGUCAUUCAGG -s CGU

import sys
import re

def main():
    import argparse
    import find_orf
    import translate

    parser = argparse.ArgumentParser(
            formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument should be treated as a path to a '
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['AUG'],
            help = ('One or more possible start codons.'))
    parser.add_argument('-x', '--stop-codons',
            type = str,
            nargs = '+', # one or more arguments
            default = ['UAA', 'UAG', 'UGA'],
            help = ('One or more possible stop codons.'))

    args = parser.parse_args()

    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    orf = find_orf.find_first_orf(sequence = sequence,
            start_codons = args.start_codons,
            stop_codons = args.stop_codons)

    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}

    translation = translate.translate_sequence(rna_sequence = orf, genetic_code = genetic_code)
    sys.stdout.write('{}\n'.format(translation))

if __name__ == '__main__':
    main()

#translate_orf.py [-h] [-p] [-s START_CODONS [START_CODONS ...]]
#                        [-x STOP_CODONS [STOP_CODONS ...]]
#                        SEQUENCE
#
#positional arguments:
#  SEQUENCE              The sequence to search for an open-reading frame. If
#                        the path flag ('-p'/'--path') is specified, then this
#                        should be a path to a file containing the sequence to
#                        be searched.
#
#optional arguments:
#  -h, --help            show this help message and exit
#  -p, --path            The sequence argument should be treated as a path to a
#                        containing the sequence to be searched. (default:
#                        False)
#  -s START_CODONS [START_CODONS ...], --start-codons START_CODONS [START_CODONS ...]
#                        One or more possible start codons. (default: ['AUG'])
#  -x STOP_CODONS [STOP_CODONS ...], --stop-codons STOP_CODONS [STOP_CODONS ...]
#                        One or more possible stop codons. (default: ['UAA',
#                        'UAG', 'UGA'])
