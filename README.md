# Make KMERS from protien fasta 
A python script to make kmers overlapping by userdefined n, additionally a PTM list is needed to tag and count the peptide combinations, see example in data folder.

A typical usage would be with following arguments.\
-fasta protein fasta file.\
-ptm An optional tab delimited file with atleast 3 columns (position, aminoacid and the ptm).\
-mer int argument specificying the length of kmer.\
-overlap int argument specifiying the overlap between subsequent kmers.\
-out a string to write the output csv.\

A example usage with PTM file would be 
```
python3 scripts/parseProt.py \
-fasta data/tauMain.fasta \
-ptm data/PTMList_DANIEL.txt \
-mer 11 \
-overlap 1 \
-out outputs/peptideList.csv
```

A example usage without PTM file would be 
```
python3 scripts/parseProt.py \
-fasta data/tauMain.fasta \
-mer 11 \
-overlap 1 \
-out outputs/NoPTMpeptideList.csv
```