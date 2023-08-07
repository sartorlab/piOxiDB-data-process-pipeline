import os,sys
import pandas as pd
import numpy as np
import subprocess
import pysam
import gzip

def read_sample_file(file):
	out = {}
	files = []
	dat = pd.read_csv(file,sep="\t",header=None)
	if np.sum(dat.count()) == dat.shape[0] * dat.shape[1]:
		with open(file) as f:
			for line in f:
				line = line.strip()
				tp = line.split("\t")
				group = tp[2]
				treat = tp[1]
				file  = tp[0]
				files.append(file)
				if group in out.keys():
					if treat in out[group].keys():
						out[group][treat].append(file)
					else:
						out[group][treat] =[]
						out[group][treat].append(file)
				else:
					out[group]={}
					if treat in out[group].keys():
						out[group][treat].append(file)
					else:
						out[group][treat] =[]
						out[group][treat].append(file)
		return(out, files,True)
	else:
		return(out,files,False)

def check_raw_files(path,samplefile,p1):
	p2 = []
	if p1 == True:
		print("checking input files ...")
		files = os.listdir(path)
		for i in files:
			if i in samplefile:
				p2.append(1)
			else:
				p2.append(0)
		if 0 in p2:
			return(False)
		else:
			return(True)
	else:
		print("Missing sample info ...")
		return(False)


def generate_output_folder(outpath,p1):
	# fastqcBefore: fastqc before trimming
	# trimming
	if p1 == True:
		print("file checked ...")
		print("mkdirs ...")
		checked_outpath = check_outpath(outpath)
		print("Output folders are creating in "+checked_outpath)
		if os.path.isdir(checked_outpath + "1_QCbeforeTrim/") is False:
			cmd1= "mkdir "+ checked_outpath + "1_QCbeforeTrim/"
			subprocess.run(cmd1,shell=True)
		if os.path.isdir(checked_outpath + "2_trimmedFastq/") is False:
			cmd2= "mkdir "+ checked_outpath + "2_trimmedFastq/"
			subprocess.run(cmd2,shell=True)
		if os.path.isdir(checked_outpath + "2_length_selection/") is False:
			cmd2s="mkdir "+ checked_outpath + "2_length_selection/"
			subprocess.run(cmd2s,shell=True)
		if os.path.isdir(checked_outpath + "3_QCafterTrim/") is False:
			cmd3= "mkdir "+ checked_outpath + "3_QCafterTrim/"
			subprocess.run(cmd3,shell=True)
		if os.path.isdir(checked_outpath + "4_Mapping/") is False:
			cmd4= "mkdir "+ checked_outpath + "4_Mapping/"
			subprocess.run(cmd4,shell=True)
		if os.path.isdir(checked_outpath + "5_Pepr/") is False:
			cmd5= "mkdir "+ checked_outpath + "5_Pepr/"
			subprocess.run(cmd5,shell=True)
		if os.path.isdir(checked_outpath + "6_piRNA/") is False:
			cmd6= "mkdir "+ checked_outpath + "6_piRNA/"
			subprocess.run(cmd6,shell=True)
		return(True)
	else:
		print("Not all fastq files existed in the data folder ...")
		return(False)

def check_outpath(outpath):
	newpath = ""
	if outpath.endswith("/"):
		newpath = outpath
	else:
		newpath = outpath+"/"
	return(newpath)

def qc_before_trim(datafd,outpath):
	checked_outpath = check_outpath(outpath)
	checked_datadf  = check_outpath(datafd)
	files = os.listdir(datafd)
	for i in files:
		cmd = "fastqc -o " + checked_outpath + "1_QCbeforeTrim/ -t 4 " + checked_datadf + i
		subprocess.run(cmd,shell=True)

def trimming_fastq_with_cutadapt(datafd,outpath):
	#
	checked_outpath = check_outpath(outpath)
	checked_datadf  = check_outpath(datafd)
	files = os.listdir(datafd)
	for i in files:
		outfile = i.rstrip("fastq.gz") + ".trimmed.fq.gz"
		cmd = "cutadapt -m 15 -u 3 -a AAAAAAAAAA -o " + checked_outpath + "2_trimmedFastq/" + outfile + " " + checked_datadf+i
		subprocess.run(cmd,shell=True)

def lengthe_selection(outpath,fd):
	checked_outpath = check_outpath(outpath)
	checked_fd      = check_outpath(checked_outpath+fd)
	files = os.listdir(checked_outpath+"2_trimmedFastq/")
	for i in files:
		if i.endswith(".trimmed.fq.gz"):
			inpfile = checked_outpath+"2_trimmedFastq/"+i
			outfile = gzip.open(checked_fd+i.rstrip("trimmed.fq.gz")+".filtered.fq.gz",'wb')
			with pysam.FastxFile(inpfile) as f:
				for entry in f:
					lens = len(entry.sequence)
					if lens <= 45 and lens >= 20:
						outfile.write((entry.name+"\n"+entry.sequence+"\n"+entry.comment+"\n"+entry.quality+"\n").encode())
			outfile.close()

def check_empty_folder(outpath,subs):
	checked_outpath = check_outpath(outpath)
	ipfd = checked_outpath + subs
	if os.path.exists(ipfd) and not os.path.isfile(ipfd):
		if not os.listdir(ipfd):
			return(False)
		else:
			return(True)
	else:
		return(False)

def qc_after_trim(outpath):
	checked_outpath = check_outpath(outpath)
	files = os.listdir(checked_outpath+"2_length_selection/")
	for i in files:
		cmd = "fastqc -o " + checked_outpath + "3_QCafterTrim/ -t 4 " + checked_outpath + "2_length_selection/" + i
		subprocess.run(cmd,shell=True)
	
def run_mapping(outpath,genome,gtf):
	checked_outpath = check_outpath(outpath)
	files = os.listdir(checked_outpath+"2_length_selection/")
	for i in files:
		tag = i.strip("fastq.gz")
		tag = tag.strip("fq.gz")
		cmd = "STAR --runThreadN 8 --outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16 --outFilterScoreMinOverLread 0  --outFilterMatchNminOverLread 0 --alignIntronMax 1 --genomeDir " + genome + " --sjdbGTFfile " + gtf + " --readFilesIn " + checked_outpath + "2_trimmedFastq/" + i + " --readFilesCommand zcat  --outSAMtype BAM SortedByCoordinate --outFileNamePrefix " + checked_outpath + "4_Mapping/STAR_sortedByCoord_" + tag+"."
		subprocess.run(cmd,shell=True)

def generate_bai_file_for_pepr(outpath):
	checked_outpath = check_outpath(outpath)
	files = os.listdir(checked_outpath+"4_Mapping/")
	for i in files:
		if i.endswith("bam"):
			cmd = "samtools index " + checked_outpath + "4_Mapping/"+i
			subprocess.run(cmd,shell=True)

def run_pepr(sample,outpath):
	for group,values in sample.items():
		untrted = values['untreated']
		treated = values['treated']
		untrtedfiles = get_bam_files(untrted,outpath)
		treatedfiles = get_bam_files(treated,outpath)
		cmd = "PePr -c "+",".join(treatedfiles) + " -i " + ",".join(untrtedfiles) + " --file-format bam --shiftsize 0 --windowsize 20 --threshold 1e-3 --peaktype sharp --num-processors 1 --name "+str(group)+" --output "+ check_outpath(outpath)+"5_Pepr/"
		subprocess.run(cmd,shell=True)


def get_bam_files(files,outpath):
	out = []
	bamfiles = os.listdir(check_outpath(outpath)+"4_Mapping/")
	for i in files:
		tag = i.strip("fastq.gz")
		tag = tag.strip("fq.gz")
		for x in bamfiles:
			if tag in x and x.endswith("bam"):
				out.append(check_outpath(outpath)+"4_Mapping/"+x)
	return(out)

def generate_piRNA_from_Peaks(outpath,fasta):
	outpath_checked = check_outpath(outpath)
	files = os.listdir(outpath_checked+"5_Pepr/")
	for i in files:
		if i.endswith("bed"):
			generate_file_with_seq(outpath_checked,fasta,i)

def get_tag(file):
	return(file.split("__")[0])

def read_fasta(fp):
	name, seq = None, []
	for line in fp:
		line = line.rstrip()
		if line.startswith(">"):
			if name: yield (name, ''.join(seq))
			name, seq = line.split()[0].strip(">"), []
		else:
			seq.append(line)
	if name: yield (name, ''.join(seq))

def generate_file_with_seq(outpath,fasta,file):
	tag = get_tag(file)
	inputfile = outpath+"5_Pepr/"+file
	outptfile = open(outpath+"6_piRNA/"+tag+".piRNA.txt",'w')
	coord = get_coord(inputfile)
	for seqid,seq in read_fasta(open(fasta,"r")):
		for chrs,v in coord.items():
			if chrs == seqid:
				for z in v:
					targetseq = seq[int(z[0])-1:int(z[1])]
					rightflnk = seq[int(z[0])-1-3:int(z[0])-1]
					leftflnak = seq[int(z[1])-1:int(z[1])-1+3]
					outline = z[2] + "\t" + rightflnk.upper() + "\t" +targetseq.upper() + "\t" + leftflnak.upper() + "\n"
					outptfile.write(outline)
	outptfile.close()


def get_coord(file):
	out = {}
	with open(file) as f:
		for line in f:
			line = line.strip()
			chrs = line.split("\t")[0]
			strt = line.split("\t")[1]
			end  = line.split("\t")[2]
			if chrs in out.keys():
				out[chrs].append([strt,end,line])
			else:
				out[chrs]=[]
				out[chrs].append([strt,end,line])
	return(out)	