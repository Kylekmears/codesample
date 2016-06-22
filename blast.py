#!/usr/bin/env python3

# In a virtualenv, run "pip install -r requirements.txt" before
# running this code

# TO DO: 
# 1) Figure out how to get sequence, description, percent ID, and E-value
#    from blast_record
# 2) Add e_value cutoff for NCBIWWW.qblast
# 3) Create sqlite database and add these to the database
# 4) Refactor so B_T_D is blast (returns xml) and xml_to_database

import argparse
import sqlite3
import os
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

def blast_to_database(in_file, out_database):
    '''Runs NCBI BLAST query against nr and stores results in an sqlite database'''

    for index, fasta_sequence in enumerate(SeqIO.parse(in_file, 'fasta')):

        print('Blasting sequence', index + 1)
        blast_record = [seq for seq in NCBIXML.parse(NCBIWWW.qblast('blastn', 'nr', fasta_sequence.seq))]
        out_file  = open('result' + str(index) + '.txt','w')
        record_copy = blast_record

        if len(blast_record) == 0:
            print('No results found')

        for blast_result in blast_record:
            print('blast_record: ', blast_record, type(blast_record))
            print('blast_result', blast_result, type(blast_result))
            out_file.write('')
            # sequence = blast_result.alignments.hsps.strand (?), description = blast_result.descriptions.title (or maybe .score), 
            # percent ID = ?, and E-value = blast_result.descriptions.e
            print([(i.num_alignments, type(i.num_alignments)) for i in blast_result.descriptions])
            # So here is what's happening: seq_record = sequence to be blasted
            # blast_record = generator that contains sequence.  sequence = blast object

def make_database(out_database):

    try:
        connection = sqlite3.connect(out_database)
        db = connection.cursor()

        db.execute("""
                   CREATE table if not exists
                   summary(title TEXT,
                           run_id INTEGER, 
                           num_hits INTEGER,
                           length INTEGER,
                           best_hit TEXT,
                           best_hit_score REAL)
                   """)

        db.execute("""
                   CREATE table if not exists
                   results(description TEXT,
                           run_id INTEGER,
                           percent_id REAL,
                           length INTEGER,
                           e_value REAL,
                           sequence TEXT)
                   """)


        connection.commit()

    except Exception as e:
        print('raised' + e + 'when trying to connect to database')

def send_to_database(out_database, blast_result):
    '''Takes in biopython Bio.Blast.Record.Blast object and
       sends important results to database'''
    
    connection = sqlite3.connect(out_database)
    db = connection.cursor()

    db.execute("""
               INSERT into summary
               (title,
                run_id,
                num_hits,
                length,
                best_hit,
                best_hit_score)
               values(VALUES)
               """)

    db.execute("""
               INSERT into results
               (description,
                run_id,
                percent_id,
                length,
                e_value,
                sequence)
               values(VALUES)
               """)

def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    
    parser.add_argument("input", help = "Path to fasta-format file containing \
                        DNA sequences", type = str)
    parser.add_argument("-o", "--output", help =
                        "Specify path to output database", type = str)
    parser.add_argument("-e", "--e_value", help =
                        "e-value cut-off. default = 0.001", type = float)

    group.add_argument("-v","--verbose", action ="store_true")
    group.add_argument("-q","--quiet", action ="store_true")
    
    args = parser.parse_args()

#    for seq_record in SeqIO.parse(args.input, 'fasta'):
#        print(seq_record.id)
#        print(repr(seq_record.seq))
#        print(len(seq_record))
#        print(seq_record.description)

    if args.output:
        out_db = args.output
    else:
        out_db = "blast_results.db"
    
    make_database(out_db)
    blast_to_database(args.input, out_db)
        
    if args.verbose:
        print('Verbose')
    elif args.quiet:
        print('shh')
        
if __name__ == "__main__":
    main()
