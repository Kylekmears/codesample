#!/usr/bin/env python3

# In a virtualenv, run "pip install -r requirements.txt" before
# running this code

import argparse
import sqlite3
import requests
import biopython

def blast(n):
    '''Web BLAST query against nr and stores results in an sqlite database'''

def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    
    parser.add_argument("input", help = "Path to fasta-format file containing \
                        DNA sequences", type = str)
    parser.add_argument("-o", "--output", help =
                        "Specify path to output database", type = str)

    group.add_argument("-v","--verbose", action ="store_true")
    group.add_argument("-q","--quiet", action ="store_true")
    
    args = parser.parse_args()

    result = fib(args.num)
    
    if args.output:
        out_db = args.output
    else:
        out_db = "blast_results.db"
        
    if args.verbose:
        print( "The " + str(args.num)+"th fib number is "+ str(result))
    elif args.quiet:
        print(result)
        
if __name__ == "__main__":
    main()
