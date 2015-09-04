# blastable
Generates tables of blast hits.


## Example usage


``python blastable.py -i multifasta.fa ../genomes/``

To display other options:

`python blastable.py -h`

## Input File Requirements

The slightly tricky bit is formatting your 
multifasta.fa headers to be compatible with BLAST.

Each fasta header must contain 4 comma separated fields. This is similar to the format used by [Seqfindr](https://github.com/mscook/seqfindr) and was designed as a follow up analysis once you have files formatted for Seqfindr.

`>field1,field2,field3,field3[category]`

* field1 = Gene id/locus tag [Optional]
* field2 = Gene name (must be unique)
* field 3 = Gene product [Optional]
* field 4 = Species [Optional]
* category = Can be any keyword used to group genes. E.g. Virulence, Resistance

Optional fields must be present even if they contain no text. E.g. `>,rpoB,,`

Before the fasta file can be used with blastable.py:

1. All spaces need to be replaced with `_` (underscore)
2. Fields that have no characters in them need to be substituted with dummy text. i.e. replace `,,` with `,dummytext,`


## Software Dependencies

1. [Command line BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/) on your $PATH
2. [Biopython](http://biopython.org/DIST/docs/install/Installation.html#sec14)
3. [Pandas](http://pandas.pydata.org/) - Python library for data analysis

If you are using Mac OSX you can use pip, a python package manager, to install `pandas`

`sudo pip install pandas`

