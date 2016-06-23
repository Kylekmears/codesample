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
from datetime import datetime

def blast_to_database(in_file, out_database):
    '''Runs NCBI BLAST query against nr and stores results in an sqlite database'''

    for index, fasta_sequence in enumerate(SeqIO.parse(in_file, 'fasta')):

        print('Blasting sequence', index + 1)
        date = str(datetime.now())
        blast_record = NCBIXML.read(NCBIWWW.qblast('blastn', 'nr', fasta_sequence.seq))

        summary_title = fasta_sequence.description
        run_id = fasta_sequence.id
        best_hit = None
        best_hit_score = 100

        for record in range(len(blast_record.descriptions)):

            top_percent = 0
            sum_len = 0
            sum_id = 0

            for hsp in blast_record.alignments[record].hsps:
                sum_id += hsp.identities
                sum_len += hsp.align_length
                hsp_percent = (hsp.identities/hsp.align_length)*100
                if hsp_percent > top_percent:
                    top_percent = hsp_percent

            description = blast_record.descriptions[record].title
            percent_id_top_HSP = top_percent
            percent_id_all_HSP = (sum_id/sum_len)*100
            e_val = blast_record.descriptions[record].e

            if best_hit_score > e_val:
                best_hit_score = e_val
                best_hit = description
            
            result = [description, run_id, percent_id_top_HSP, percent_id_all_HSP, e_val, date]
            result_to_database(out_database, result)

        num_hits = len(blast_record.descriptions)
        if num_hits == 0:
            print('No results found')
            
        summary = [summary_title, run_id, num_hits, best_hit, best_hit_score, date]
        summary_to_database(out_database, summary)

def make_database(out_database):

    try:
        connection = sqlite3.connect(out_database)
        db = connection.cursor()

        db.execute("""
                   CREATE table if not exists
                   summary(title TEXT,
                           run_id TEXT, 
                           num_hits INTEGER,
                           top_hit TEXT,
                           top_hit_e_value REAL,
                           date TEXT)
                   """)

        db.execute("""
                   CREATE table if not exists
                   results(description TEXT,
                           run_id TEXT,
                           percent_id_top_HSP REAL,
                           percent_id_all_HSP REAL,
                           e_value REAL,
                           date TEXT)
                   """)


        connection.commit()

    except Exception as e:
        print('raised' + e + 'when trying to connect to database')


def summary_to_database(out_database, summary):
    '''Takes in summary list of results and sends data to database
        overview = [title, run id, num hits, best hit, best hit score]'''

    connection = sqlite3.connect(out_database)
    db = connection.cursor()

    db.execute("""
               INSERT into summary
               (title, run_id, num_hits, top_hit, top_hit_e_value, date)
               values(?, ?, ?, ?, ?, ?)
               """, summary
              )

    connection.commit()

def result_to_database(out_database, result):
    '''takes in results list and adds data to results table
       result = [description, run id, percent id top HSP, percent id all HSP, e val]'''

    connection = sqlite3.connect(out_database)
    db = connection.cursor()
    db.execute("""
               INSERT into results
               (description, run_id, percent_id_top_HSP, 
                percent_id_all_HSP, e_value, date)
               values(?, ?, ?, ?, ?, ?)
               """, result
              )

    connection.commit()

def main():

    parser = argparse.ArgumentParser()
    
    parser.add_argument("input", help = "Path to fasta-format file containing \
                        DNA sequences", type = str)
    parser.add_argument("-o", "--output", help =
                        "Specify path to output database (directories \
                         in path must already exist)", type = str)

    args = parser.parse_args()

    if args.output:
        out_db = args.output
    else:
        out_db = "blast_results.db"
    
    make_database(out_db)
    blast_to_database(args.input, out_db)
    print('Complete. Results sent to', out_db)
        
if __name__ == "__main__":
    main()
