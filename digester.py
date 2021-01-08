#!/usr/bin/env python3

from Bio.Seq import Seq
from Bio import SeqIO, Restriction
from Bio.Restriction import *
import pandas as pd
import argparse
import re

def parse_seq(fasta):
    """Parse fasta format file with DNA sequences. """
    fasta = open(fasta)
    sequences = []
    for record in SeqIO.parse(fasta, "fasta"):
        dna = {'id' : record.id, 'seq' : record.seq}
        sequences.append(dna)
    return sequences

def enzymes_list(enzymes):
    """Imports list of restriction enzymes from the file."""
    with open(enzymes) as f:
        enzyme_list = f.readlines()
        for line in enzyme_list:
            enzymes = line.split(",")
    rest_batch = RestrictionBatch(enzymes)
    return rest_batch

def fragment_lenght(enzyme, sequence):
    """Calculates lenght of DNA fragments after restriction digest."""
    digest = enzyme.catalyse(sequence, linear=args.circular)
    lenght = [len(d) for d in digest]
    return lenght

def digest():
    """Do the restriction digest"""
    enz = enzymes_list(args.enzymes)
    sequences = parse_seq(args.input)
    restriction_digest = []
    for s in sequences:
        for enzyme in enz:
            search_results = {
                'id': s['id'], 
                'enzyme': str(enzyme), 
                'location': str(enzyme.search(s['seq'], linear=args.circular)),
                'fragments': fragment_lenght(enzyme, s['seq']) # function fragment_lenght() above
                }
            restriction_digest.append(search_results)
    return restriction_digest

def results():
    """Saves results of restriction digest to a txt and csv file."""
    restriction_digest = digest()
    filename = open(args.output, "w")
    for r in restriction_digest:
        if len(r['location']) > 2:
            filename.write("{} digested with {}".format(r['id'], r['enzyme']))
            filename.write("\nRestriction sites location: {}".format(r['location']))
            filename.write("\nFragments: {}\n\n".format(str(r['fragments'])))
        else:
            filename.write("{} no restriction sites for {}\n \n".format(r['id'], r['enzyme']))
    filename.close()
    df = pd.DataFrame(restriction_digest)
    df.columns = ['enzyme','fragments','sequence name', 'location']
    df.to_csv(f'{args.output}.csv')

def verbose():
    seq_num = parse_seq(args.input)
    enzymes = enzymes_list(args.enzymes)
    enz = re.split('\+', str(enzymes)) 
    print("\nDigested {} sequences with: \n".format(len(seq_num)))
    for e in enzymes:
        print(e)
        #print(("{} " * len(enz)).format(*enz))
    if args.circular == False:
        print("\nCircular sequence option was selected.")
    else:
        print("Linear sequence option was selected.")

def make_parser():
    parser = argparse.ArgumentParser(prog='Digester.py', description="This script performs restriction digest of DNA.")
    parser.add_argument('-i', '--input', help="Fasta file storing sequences you would like to digest", type=str, required=True)
    parser.add_argument('-o','--output', help="Output file name", type=str, required=True)
    parser.add_argument('-e', '--enzymes',help="Restriction enzymes list", type=str, required=True)
    parser.add_argument('-v', '--verbose', action='store_true', help="Print number of sequences and restriction enzymes used.")
    parser.add_argument('-c', '--circular', action='store_false', help="Circular sequence, e.g. plasmid") # controls option in fragment_lenght() function
    return parser

if __name__ == '__main__':
    parser = make_parser()
    args = parser.parse_args()

    if args.verbose:
        verbose()
    
    results()

