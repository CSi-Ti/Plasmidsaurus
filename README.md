
# Plasmidsaurus QC

### System requirement
```text
Linux, Windows (via WSL), [MacOS not tested]
```

### Environment installation

```text
# install conda

# install conda environment
conda env create -f ./environment.yml
```

### File setting
```text
- modify the last row in run.py to specify your fastq folder and reference file.
- the reference file should have two columns: TCR_alpha_full, TCR_beta_full
```

### Execution
```text
conda run -n plasmidsaurus --no-capture-output python ./run.py 
```

### Output
```text
Result will appear in ./data/result.csv
The columns are:
- File_name: file name of fastq
- Total_reads: total read number in fastq
- Kozak_reads: filtered read number with kozak pattern
- Type: TCR-alpha or TCR-beta
- Type_reads: filtered read number with TCR alpha/beta pattern
- Lib_seq: TCR nucleotide sequence from reference file
- Cons_seq: consensus sequence derived from fastq
- Match: if Lib_seq equals to Cons_seq
- N_match: number of reads which can perfectly match to Cons_seq
- Percentage: N_match / Type_reads
```