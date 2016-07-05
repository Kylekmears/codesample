Module requirements can be found in requirements.txt and installed with the following command:

pip3 install -r requirements.txt

_____________________________________
usage: blast.py [-h] [-o OUTPUT] input

positional arguments:
  input                 Path to fasta-format file containing DNA sequences.

optional arguments:
  -h, --help            show this help message and exit
  
  -o OUTPUT, --output OUTPUT
  
                        Specify path to output database (directories in path
                        must already exist)

_________________________________________
This code results in a database with two tables:
summary, which contains title, run id, number of hits, best hit, and best hit e-value.  
And results, containing description, run id, percent id for top HSP, percent id for all HSPs, and e val

The default database name is blast_database.db
