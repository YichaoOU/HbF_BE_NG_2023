# Summary

The bam files were generated using the RNA-seq variant calling pipeline here: https://hemtools.readthedocs.io/en/latest/content/NGS_pipelines/rna_seq_variant_call.html. Log files can be found in: `rna_seq_variant_call_yli11_2023-01-09`

Then the bam files were analyzed using `REDItools`: https://github.com/BioinfoUNIBA/REDItools. Specifically, the first step (`run.lsf`) is:

```
python ~/HemTools/share/script/REDItoolDnaRna.py -t 20 -i ${COL1} -d -D -f /home/yli11/Data/Human/hg19/fasta/Homo_sapiens.add_chr.GRCh37.dna.primary_assembly.fa -G ~/Data/Human/hg19/annotations/gencode.v30lift37.annotation.st.gtf.gz -o ${COL1}_REDItool_results
```

And the second step (`run2.lsf`) is to:

```

python ~/HemTools/share/script/tabulate_reads_AtoI-stranded.py ${COL1} 

```

The deamination frequency is then generated in: `analysis.ipynb`. 

REF: https://www.nature.com/articles/s41587-020-0453-z


