import os,sys
from os import path
import argparse
from functions import *
import subprocess

# this script is to process Sodium periodate sequencing data for pioxiDB
# 1. fastqc before trimming (including multiqc)
# 2. trimming
######length selection 
# 3. fastqc after  trimming (including multiqc)
# 4. mapping (genome should be indexed)
# 5. Peak calling
# 6. change format for data dump

def main(argv):
	# subprocess.run("ml python",shell=True)
	# subprocess.run("ml star",shell=True)
	parser = argparse.ArgumentParser(description="This code is used to process Sodium periodate RNA sequencing data and change the format for the database.\n Example:\n python pioxidata.process.py -f path/to/fastq/files -s sampleinfoFile -o path/to/output/folder -g path/to/genome/index -m False -t GTFFile -G genomeFastaFile")
	parser.add_argument("-f","--folderPath",help="This is the folder that contains the raw fastq files")
	parser.add_argument("-s","--sampleInfo",help="This is the sample info file, two columns without header, 1st: sample ID; 2nd: treatment(i.e., treated or untreated); 3rd: group; tsv format. See example file sample.info.txt for details.")
	parser.add_argument("-o","--outputPath",help="Output folder. if not set, output will be generate in the current folder.")
	parser.add_argument("-g","--genomePath",help="Genome index location, need index the genome and provide the GTF file.")
	parser.add_argument("-m","--multiqc",help="run multiqc or not; default to False")
	parser.add_argument("-t","--gtf",help="GTF file for mapping")
	parser.add_argument("-G","--genomefasta",help="fasta file of reference genome for extracting piRNA sequences")
	
	args=parser.parse_args()

	if args.folderPath == None:
		print("No data folder provided!")
		exit()
	else:
		datafd=args.folderPath

	if args.sampleInfo == None:
		print("No sample info provided! Check help for more info!")
		exit()
	else:
		sample=args.sampleInfo

	if args.outputPath == None:
		print("Generating results in current folder ...")
		outpath= os.path.dirname(os.path.abspath(__file__)) + "/"
	else:
		outpath=args.outputPath

	if args.genomePath == None:
		print("No genome index path provided!")
		exit()
	else:
		genome=args.genomePath

	if args.gtf == None:
		print("No GTF provided!")
		exit()
	else:
		gtf=args.gtf

	if args.genomefasta == None:
		print("No genome fasta provided!")
		exit()
	else:
		genomefa = args.genomefasta

	if args.multiqc == None:
		multi = False
	elif args.multiqc == "T" or args.multiqc == "True" or args.multiqc == "true":
		multi = True
	elif args.multiqc == "F" or args.multiqc == "False" or args.multiqc == "False":
		multi = False
	else:
		multi = False

	# check fastq files
	sampledict,files,p1 = read_sample_file(sample)
	p2=check_raw_files(datafd,files,p1)
	# generate output folder
	p3=generate_output_folder(outpath,p2)
	# start running the pipeline
	# 1 qc before trimming
	p4=check_empty_folder(outpath,"1_QCbeforeTrim/")
	if p4 is True:
		print("Qc is done ...")
	else:	
		print("QC is running ...")
		qc_before_trim(datafd,outpath)
		print("Qc is done ...")
	# 2 trimming
	p5=check_empty_folder(outpath,"2_trimmedFastq/")
	if p5 is True:
		print("trimming is done ...")
	else:
		print("start trimming ...")
		trimming_fastq_with_cutadapt(datafd,outpath) # this function might need to be changed for different sequencing data
		print("trimming is done ...")
	p5s=check_empty_folder(outpath,"2_length_selection/")
	if p5s is True:
		print("length selection is done ...")
	else:
		print("selecting reads ...")
		lengthe_selection(outpath,"2_length_selection/")
		print("length selection is done ...")
	# 3 QC after trimming
	p6 = check_empty_folder(outpath,"3_QCafterTrim/")
	if p6 is True:
		print("Qc after trimming is done ...")
	else:
		print("QC is running ...")
		qc_after_trim(outpath)
		print("Qc after trimming is done ...")
	# 4 mapping: need indexed genome and corresponding GTF files 
	p7 = check_empty_folder(outpath,"4_Mapping/")
	if p7 is True:
		print("Mapping is done ...")
	else:
		print("Mapping is running ...")
		run_mapping(outpath,genome,gtf)
		print("Mapping is done ...")
	# 5 PePr: peak calling
	p8 = check_empty_folder(outpath,"5_Pepr/")
	if p8 is True:
		print("Pepr is done ...")
	else:
		print("PePr is running ...")
		generate_bai_file_for_pepr(outpath)
		run_pepr(sampledict,outpath)
		print("Pepr is done ...")

	p9 = check_empty_folder(outpath,"6_piRNA/")
	if p9 is True:
		print("piRNA generated ...")
	else:
		print("generating piRNA ...")
		generate_piRNA_from_Peaks(outpath,genomefa)
		print("piRNA generated ...")

if __name__ =="__main__":
	main(sys.argv[1:])