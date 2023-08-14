# piOxiDB-data-process-pipeline

This repository contains the data process pipeline for identifying piRNAs from sodium periodate treatment small RNA data.
The piOxiDB can be found [here](https://pioxidb.dcmb.med.umich.edu/).

This pipeline is written in Python 3.X.

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://multiqc.info/), [Cutadapt](https://cutadapt.readthedocs.io/en/stable/), [STAR](https://github.com/alexdobin/STAR), and [PePr](https://github.com/shawnzhangyx/PePr) are needed to run this Python code. 

See more details, run:

```
python pioxidata.process.py -h
```