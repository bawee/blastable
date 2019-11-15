# blastable
Generates tables of blast hits.


*Note:* The default threshold is 60% average sequence identity of the alignment divided by the length of the query. 

`thresh = ((pident*HSPlength)/qlen)`

## Example usage


``python blastable.py -i multifasta.fa ../genomes/``

To display other options:

`python blastable.py -h`

## Input File Requirements

The slightly tricky bit is formatting your 
multifasta.fa headers to be compatible with BLAST.

Each fasta header must contain 2 comma separated fields.

`>field1,field2`

* field1 = Gene name (must be unique)
* field2 = Description (free text, no commas please)

Before the query fasta file can be used with blastable.py:

1. All spaces should to be replaced with `_` (underscore)

## Software Dependencies

1. [Command line BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/) on your $PATH
2. [Biopython](http://biopython.org/DIST/docs/install/Installation.html#sec14)
3. [Pandas](http://pandas.pydata.org/) - Python library for data analysis

You can use pip, a python package manager, to install `pandas`

`sudo pip install pandas`

## Version History


* v0.1 - Initial version
* v0.2 - Sorts blast results based on the input sequence. (28 Sep 2015)
* v0.3 - Handles empty columns. I.e. Queries with no hits in any of the genomes (28 Sep 2015).
* v0.4 - Prints BLAST results into a separate folder within the working directory. (28 Sep 2015)
* v0.5 - Also performs blastp. Removed seqfindr format requirement (17 May 2018)

