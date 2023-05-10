# Summary

The CUT & RUN results were generated using `HemTools cut_run -f fastq.tsv -p peakcall.tsv`. Log files can be found in the `log_files` folder.

The pipeline is written in a Python script `HemTools`, which calls the sub-command `cut_run.py`. 

Dependencies:

- fastqc

- bwa

- samtools

- deeptools

- macs2

- bedtools