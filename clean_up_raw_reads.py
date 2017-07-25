#	Copyright (C) 2017 Mozes Blom
#
#	This is a script to clean-up PAIRED-END shotgun sequence data. It is designed to deal with sequencing libraries from both fresh and/or museum tissues.
#
#	Please direct all questions to: mozes.blom@gmail.com
#
#	NOTE: if running on cluster, occasionally one can run into a EoF error when using the gzip python module.

#########################
## Required modules
#########################
try:
	import os
	import csv
	import subprocess
	import sys
	import datetime
	import gzip
except ImportError:
	sys.exit("One of the required modules can't be found...")

#########################
## Input
#########################
# Run specific folder name, where the results will be stored [char] 
results_path = '/Path/to/existing/folder/'
# Folder where all the raw sequence reads are stored [char] 
reads_path = '/Path/to/folder/with/reads/' 
## Path to list with individuals, list should be a Tab-delim text file that contains a header with 'lib' and 'type'. Type = 'SKIN' or 'FRESH' [CHAR]
## Expects libs to be in format: indivA_Lxxx. In this case the read file name should be: indivA_Lxxx.fastq.gz
## If 'SKIN', then the libraries are assumed to have been constructed based on DNA from museum skins and thus most F and R have been merged
## If 'FRESH', then the libraries are assumed to have been constructed based on DNA from fresh tissues and thus read pairs are still retained rather than merged
lib_list = '/Path/to/lib_list.txt'

	
#########################
# Which Analysis to run? 1 = yes, 0 = no
#########################

# Clean up data
complete = 0       			# [int, 0 or 1] Run complete analysis?

# ----------------------- #

# (Re)run subset of analyses; note expects the same folder structure, as if running complete analyses! Mate pair data will not be merged since there should be no overlap
pre_fastqc = 0				# [int, 0 or 1] Run FastQC on the raw reads
dedup = 0					# [int, 0 or 1] Run deduplication using super deduper
trim_adapt = 0				# [int, 0 or 1] Remove sequencing adapters with Trimmomatic, no quality control!!
merge = 0					# [int, 0 or 1] Merge reads with x overlap
trim_qual = 0				# [int, 0 or 1] Remove low quality fragments
remove_complex = 1			# [int, 0 or 1] Remove low complexity reads

# ----------------------- #

# Optional
limit_storage = 1			# [int, 0 or 1] Only store final cleaned up reads and remove all the intermediate folders, to save disk space!


#########################
## Path to dependables
#########################

fastqc = '/Path/to/bin/FastQC/fastqc'										# [char] Path to compiled FastQC executable
superdedup = '/Path/to/bin/super_deduper-1.4/super_deduper'					# [char] Path to compiled Superdeduper executable
trimmomatic = '/Path/to/bin/Trimmomatic-0.36/trimmomatic-0.36.jar'			# [char] Path to compiled Trimmomatic executable
pear = '/Path/to/bin/pear-0.9.10-bin-64/pear-0.9.10-bin-64'					# [char] Path to compiled PEAR executable
pigz = '/Path/to/bin/pigz-2.3.4/pigz'										# [char] Path to compiled PIGZ executable

#########################
## Parameter settings
#########################

threads = 8						# [int] Number of available threads to use simultaneously, i.e. for FastQC, pigz, etc.
threads_max = 8					# [int] Number of available threads to use for computational bottle necks (i.e. PEAR)
ram = 48						# [int] Max. memory available, in Gb
min_read_len = 30				# [int] Minimum length of reads (or merged reads) to be retained
min_overlap = 20				# [int] Minimum overlap between pair to be considered for merging
merge_hit_threshold = 0.0001	# [int] Probability threshold for reads to be merged yes/no, see PEAR documentation (Note; default is 0.01, but I noticed a relatively large proportion of false positives: reads that were merged, but should not have been merged)
cut_off_complexity = 0.5		# [int, min 0 - max 1] Max. relative frequency of A/T/C/G/N in a read to be considered low complexity. I.e. 0.5 means all reads that contain more than half A's are removed
adapter_trim_file = '/home/mobl/bin/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa' 		# [char] Path to adapter file to be used by Trimmomatic

#########################
## Functions
#########################

## You can make modifications to the functions, i.e. change specific flags like memory or threads, just make sure the required number of arguments and order of arguments stay the same.

# Run FastQC
def runFastQC_raw(fastqc, threads, out_dir, in_file_path_r1, in_file_path_r2):
	subprocess.call("%s -t %s --outdir %s %s %s" % (fastqc, threads, out_dir, in_file_path_r1, in_file_path_r2), shell=True)

# Run Superdeduper for deduplication
def runDedup(superdedup, in_file_path_r1, in_file_path_r2, prefix):
	subprocess.call("%s -1 %s -2 %s -p %s" % (superdedup, in_file_path_r1, in_file_path_r2, prefix), shell=True)

# Run Trimmomatic to trim adapters
def runTrimmo_adapt_PE(trimmomatic, ram, threads, reads_1_in, reads_2_in, reads_1P_out, reads_1U_out, reads_2P_out, reads2U_out, adapter_file):
	subprocess.call("java -Xmx%sG -jar %s PE -threads %s %s %s %s %s %s %s ILLUMINACLIP:%s:2:30:10:8:TRUE" % (str(ram),trimmomatic, threads, reads_1_in, reads_2_in, reads_1P_out, reads_1U_out, reads_2P_out, reads2U_out, adapter_file), shell=True)

# Run Trimmomatic to trim quality
def runTrimmo_qual_PE(trimmomatic, ram, threads, reads_1_in, reads_2_in, reads_1P_out, reads_1U_out, reads_2P_out, reads2U_out, min_read_len):
	subprocess.call("java -Xmx%sG -jar %s PE -threads %s %s %s %s %s %s %s LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:%s" % (str(ram), trimmomatic, threads, reads_1_in, reads_2_in, reads_1P_out, reads_1U_out, reads_2P_out, reads2U_out, min_read_len), shell=True)

# Run Trimmomatic to trim quality
def runTrimmo_qual_SE(trimmomatic, ram, threads, reads_in, reads_out, min_read_len):
	subprocess.call("java -Xmx%sG -jar %s SE -threads %s %s %s LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:%s" % (str(ram), trimmomatic, threads, reads_in, reads_out, min_read_len), shell=True)

# Run PEAR to merge overlapping read pairs
def runPEAR(pear, reads_1_in, reads_2_in, base_name, pvalue, min_overlap, min_read_len, threads):
	subprocess.call("%s -f %s -r %s -o %s -p %s -v %s -n %s -k -j %s" % (pear, reads_1_in, reads_2_in, base_name, pvalue, min_overlap, min_read_len, threads), shell=True)

# Identify Low Complexity reads from .fq.gz file; low complexity reads have more than x percent (cut_off) identical bases
def flagLowComplexity(fq_gz_file, cut_off):
	read_set = set()
	with gzip.open(fq_gz_file, 'rb') as fq_reads:
		while 1:
			read_name = fq_reads.readline().strip().split(' ')[0]
			if not read_name: break
			seq = fq_reads.readline().strip()
			fq_reads.readline()
			fq_reads.readline()
			seq_len = float(len(seq))
			for base in ['A','C','T','G','N']:
				if float(seq.count(base))/seq_len >= cut_off:
					read_set.add(read_name)
	return read_set

def gzip_fn(pigz, threads, path_to_fn):
	subprocess.call("%s -p %s -v -r %s" % (pigz, threads, path_to_fn), shell=True)

def getReadStats(fq_gz_file):
	with gzip.open(fq_gz_file, 'rb') as fq_reads:
		tot_len = 0
		tot_num = 0
		while 1:
			read_name = fq_reads.readline().strip().split(' ')[0]
			if not read_name: break
			seq = fq_reads.readline().strip()
			tot_len += float(len(seq))
			tot_num += 1
			fq_reads.readline()
			fq_reads.readline()
		avg_len = float(tot_len)/float(tot_num)
		return tot_num, avg_len


#########################
## Automated Analysis
#########################


# Print the date + time, when analysis started
print('Clean up started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
sys.stdout.flush()

# Create the output folder in case it doesn't exist yet
if (complete == 1) or (pre_fastqc == 1):
	if os.path.isdir(results_path):
		pass
	else:
		subprocess.call("mkdir '%s'" % (results_path), shell=True)
else:
	pass


# Create a list of all the libraries to be run
# Expects libs to be in format: indivA_Lxxx. In this case the read file name would be: indivA_Lxxx.fastq.gz
libs = []
lib_type = {}
try:
	with open(lib_list, 'r') as fn:
		reader = csv.DictReader(fn, delimiter = '\t')
		for line in reader:
			lib_type[line['lib']] = line['type']
			libs.append(line['lib'])
except IOError:
	sys.exit("Can't open the file with the library names to run, did you specify the path to the file correctly; see READS_LIST in config script?")
except:
	sys.exit("Unexpected error...")


# Run FastQC on the raw sequence reads
outDir_1 = os.path.join(results_path, "1.Pre_clean_FastQC")
if (complete == 1) or (pre_fastqc == 1):
	if os.path.isdir(outDir_1):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_1), shell=True)
	lib_stats_out_F = os.path.join(outDir_1, 'pre_clean_F.csv')
	lib_stats_out_R = os.path.join(outDir_1, 'pre_clean_R.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	# generate a variable that keeps track whether the libs are 
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_1_indiv = os.path.join(outDir_1, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_1_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_1_indiv_lib = os.path.join(outDir_1_indiv, lib)
		if os.path.isdir(outDir_1_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_1_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_1_indiv_lib), shell=True)
		os.chdir(outDir_1_indiv_lib)
		lib_1_path_gz = os.path.join(reads_path, (lib + '_R1.fastq.gz'))
		lib_2_path_gz = os.path.join(reads_path, (lib + '_R2.fastq.gz'))
		try:
			runFastQC_raw(fastqc, threads, outDir_1_indiv_lib, lib_1_path_gz, lib_2_path_gz)
		except:
			sys.exit("could not execute FastQC, check input path of files and/or path to FastQC executable")
		x, y = getReadStats(lib_1_path_gz)
		writer_f.writerow((lib, x, y))
		x, y = getReadStats(lib_2_path_gz)
		writer_r.writerow((lib, x, y))
	fh_f.close()
	fh_r.close()
else:
	pass


# Trim away adapter sequences
outDir_2 = os.path.join(results_path, "2.deduplication")
if (complete == 1) or (dedup == 1):
	# Print the date + time, when trimming started
	print('Deduplication of reads started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	sys.stdout.flush()
	if os.path.isdir(outDir_2):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_2), shell=True)
	lib_stats_out_F = os.path.join(outDir_2, 'post_dedup_F.csv')
	lib_stats_out_R = os.path.join(outDir_2, 'post_dedup_R.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_2_indiv = os.path.join(outDir_2, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_2_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_2_indiv_lib = os.path.join(outDir_2_indiv, lib)
		if os.path.isdir(outDir_2_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_2_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_2_indiv_lib), shell=True)
		os.chdir(outDir_2_indiv_lib)
		lib_1_path_gz = os.path.join(reads_path, (lib + '_R1.fastq.gz'))
		lib_2_path_gz = os.path.join(reads_path, (lib + '_R2.fastq.gz'))
		# Start deduplication!!
		runDedup(superdedup, lib_1_path_gz, lib_2_path_gz, lib)
		# Compress with Pigz
		lib_1_path_dedup = os.path.join(outDir_2_indiv_lib, (lib + '_nodup_PE1.fastq'))
		lib_2_path_dedup = os.path.join(outDir_2_indiv_lib, (lib + '_nodup_PE2.fastq'))
		gzip_fn(pigz, threads, lib_1_path_dedup)
		gzip_fn(pigz, threads, lib_2_path_dedup)
		x, y = getReadStats(os.path.join(outDir_2_indiv_lib, (lib + '_nodup_PE1.fastq.gz')))
		writer_f.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_2_indiv_lib, (lib + '_nodup_PE2.fastq.gz')))
		writer_r.writerow((lib, x, y))
	fh_f.close()
	fh_r.close()


# Trim away adapter sequences
outDir_3 = os.path.join(results_path, "3.adapter_trimmed")
if (complete == 1) or (trim_adapt == 1):
	# Print the date + time, when trimming started
	print('Removal of adapter sequences using Trimmomatic started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	sys.stdout.flush()
	if os.path.isdir(outDir_3):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_3), shell=True)
	lib_stats_out_F = os.path.join(outDir_3, 'post_adapt_trim_F.csv')
	lib_stats_out_R = os.path.join(outDir_3, 'post_adapt_trim_R.csv')
	lib_stats_out_U = os.path.join(outDir_3, 'post_adapt_trim_U.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	fh_u = open(lib_stats_out_U, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_u = csv.writer(fh_u)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_u.writerow(('Libary', 'Number of reads', 'Average length per read'))
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_3_indiv = os.path.join(outDir_3, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_3_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_3_indiv_lib = os.path.join(outDir_3_indiv, lib)
		if os.path.isdir(outDir_3_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_3_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_3_indiv_lib), shell=True)
		os.chdir(outDir_3_indiv_lib)
		lib_1_path_gz = os.path.join(outDir_2, indiv, lib, (lib + '_nodup_PE1.fastq.gz'))
		lib_2_path_gz = os.path.join(outDir_2, indiv, lib, (lib + '_nodup_PE2.fastq.gz'))
		runTrimmo_adapt_PE(trimmomatic, ram, threads, lib_1_path_gz, lib_2_path_gz, (lib + '_1_p.fastq.gz'), (lib + '_1_u.fastq.gz'), (lib + '_2_p.fastq.gz'), (lib + '_2_u.fastq.gz'), adapter_trim_file)
		x, y = getReadStats((lib + '_1_p.fastq.gz'))
		writer_f.writerow((lib, x, y))
		x, y = getReadStats((lib + '_2_p.fastq.gz'))
		writer_r.writerow((lib, x, y))
		a, b = getReadStats((lib + '_1_u.fastq.gz'))
		c, d = getReadStats((lib + '_2_u.fastq.gz'))
		x = (float(a) + float(c))
		y = (float(b) + float(d))/2
		writer_u.writerow((lib, x, y))
		# Remove the dedup reads if limit storage == 1
		if limit_storage == 1:
			print('Start removing dedup reads lib %s' % (lib))
			sys.stdout.flush()
			subprocess.call("rm %s" % (lib_1_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_2_path_gz), shell=True)
			print('Finalized removing dedup reads lib %s' % (lib))
			sys.stdout.flush()
		else:
			pass
	fh_f.close()
	fh_r.close()
	fh_u.close()

# Merge overlapping read pairs
outDir_4 = os.path.join(results_path, "4.merged_reads")
if (complete == 1) or (merge == 1):
	# Print the date + time, when merging started
	print('Merging reads started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	sys.stdout.flush()
	if os.path.isdir(outDir_4):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_4), shell=True)
	lib_stats_out_F = os.path.join(outDir_4, 'post_merge_F.csv')
	lib_stats_out_R = os.path.join(outDir_4, 'post_merge_R.csv')
	lib_stats_out_U = os.path.join(outDir_4, 'post_merge_U.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	fh_u = open(lib_stats_out_U, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_u = csv.writer(fh_u)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_u.writerow(('Libary', 'Number of reads', 'Average length per read'))
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_4_indiv = os.path.join(outDir_4, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_4_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_4_indiv_lib = os.path.join(outDir_4_indiv, lib)
		if os.path.isdir(outDir_4_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_4_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_4_indiv_lib), shell=True)
		os.chdir(outDir_4_indiv_lib)
		# Merge PE reads if needed using PEAR
		lib_R1_P_path_gz = os.path.join(outDir_3, indiv, lib, (lib + '_1_p.fastq.gz'))
		lib_R2_P_path_gz = os.path.join(outDir_3, indiv, lib, (lib + '_2_p.fastq.gz'))
		lib_R1_U_path_gz = os.path.join(outDir_3, indiv, lib, (lib + '_1_u.fastq.gz'))
		lib_R2_U_path_gz = os.path.join(outDir_3, indiv, lib, (lib + '_2_u.fastq.gz'))		
		runPEAR(pear, lib_R1_P_path_gz, lib_R2_P_path_gz, lib, merge_hit_threshold, min_overlap, min_read_len, threads_max)
		# Compress PEAR output files, rename accordingly and merge the unpaired reads
		gzip_fn(pigz, threads, os.path.join(outDir_4_indiv_lib, (lib + '.unassembled.forward.fastq')))
		gzip_fn(pigz, threads, os.path.join(outDir_4_indiv_lib, (lib + '.unassembled.reverse.fastq')))
		gzip_fn(pigz, threads, os.path.join(outDir_4_indiv_lib, (lib + '.assembled.fastq')))
		subprocess.call("mv %s %s" % (os.path.join(outDir_4_indiv_lib, (lib + '.unassembled.forward.fastq.gz')), os.path.join(outDir_4_indiv_lib, (lib + '_1_p.fastq.gz'))), shell=True)
		subprocess.call("mv %s %s" % (os.path.join(outDir_4_indiv_lib, (lib + '.unassembled.reverse.fastq.gz')), os.path.join(outDir_4_indiv_lib, (lib + '_2_p.fastq.gz'))), shell=True)
		subprocess.call("cat %s %s %s > %s" % (os.path.join(outDir_4_indiv_lib, (lib + '.assembled.fastq.gz')), lib_R1_U_path_gz, lib_R2_U_path_gz, os.path.join(outDir_4_indiv_lib, (lib + '_u.fastq.gz'))), shell=True)
		x, y = getReadStats(os.path.join(outDir_4_indiv_lib, (lib + '_1_p.fastq.gz')))
		writer_f.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_4_indiv_lib, (lib + '_2_p.fastq.gz')))
		writer_r.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_4_indiv_lib, (lib + '_u.fastq.gz')))
		writer_u.writerow((lib, x, y))
		# Remove the trimmomatic reads if limit storage == 1
		if limit_storage == 1:
			print('Start removing reads lib %s' % (lib))
			sys.stdout.flush()
			subprocess.call("rm %s" % (lib_R1_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_R2_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_R1_U_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_R2_U_path_gz), shell=True)
			subprocess.call("rm %s" % (os.path.join(outDir_4_indiv_lib, (lib + '.assembled.fastq.gz'))), shell=True)
			subprocess.call("rm %s" % (os.path.join(outDir_4_indiv_lib, (lib + '.discarded.fastq'))), shell=True)
			print('Finalized removing reads lib %s' % (lib))
			sys.stdout.flush()
		else:
			pass
	fh_f.close()
	fh_r.close()
	fh_u.close()

# Trim away low-quality reads
outDir_5 = os.path.join(results_path, "5.quality_trimmed")
if (complete == 1) or (trim_qual == 1):
	# Print the date + time, when merging started
	print('Quality trimming of reads started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	sys.stdout.flush()
	if os.path.isdir(outDir_5):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_5), shell=True)
	lib_stats_out_F = os.path.join(outDir_5, 'post_qual_trim_F.csv')
	lib_stats_out_R = os.path.join(outDir_5, 'post_qual_trim_R.csv')
	lib_stats_out_U = os.path.join(outDir_5, 'post_qual_trim_U.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	fh_u = open(lib_stats_out_U, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_u = csv.writer(fh_u)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_u.writerow(('Libary', 'Number of reads', 'Average length per read'))
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_5_indiv = os.path.join(outDir_5, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_5_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_5_indiv_lib = os.path.join(outDir_5_indiv, lib)
		if os.path.isdir(outDir_5_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_5_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_5_indiv_lib), shell=True)
		os.chdir(outDir_5_indiv_lib)
		lib_R1_P_path_gz = os.path.join(outDir_4, indiv, lib, (lib + '_1_p.fastq.gz'))
		lib_R2_P_path_gz = os.path.join(outDir_4, indiv, lib, (lib + '_2_p.fastq.gz'))
		runTrimmo_qual_PE(trimmomatic, ram, threads, lib_R1_P_path_gz, lib_R2_P_path_gz, (lib + '_1_p.fastq.gz'), (lib + '_1_u.fastq.gz'), (lib + '_2_p.fastq.gz'), (lib + '_2_u.fastq.gz'), min_read_len)
		lib_U_path_gz = os.path.join(outDir_4, indiv, lib, (lib + '_u.fastq.gz'))
		runTrimmo_qual_SE(trimmomatic, ram, threads, lib_U_path_gz, (lib + '_unpaired.fastq.gz'), min_read_len)
		# merge the unpaired reads
		subprocess.call("cat %s %s %s > %s" % ((lib + '_1_u.fastq.gz'), (lib + '_2_u.fastq.gz'), (lib + '_unpaired.fastq.gz'), (lib + '_u.fastq.gz')), shell=True)
		x, y = getReadStats(os.path.join(outDir_5_indiv_lib, (lib + '_1_p.fastq.gz')))
		writer_f.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_5_indiv_lib, (lib + '_2_p.fastq.gz')))
		writer_r.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_5_indiv_lib, (lib + '_u.fastq.gz')))
		writer_u.writerow((lib, x, y))
		if limit_storage == 1:
			print('Start removing reads lib %s' % (lib))
			sys.stdout.flush()
			subprocess.call("rm %s" % (lib_R1_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_R2_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_U_path_gz), shell=True)
			subprocess.call("rm %s" % (os.path.join(outDir_5_indiv_lib, (lib + '_1_u.fastq.gz'))), shell=True)
			subprocess.call("rm %s" % (os.path.join(outDir_5_indiv_lib, (lib + '_2_u.fastq.gz'))), shell=True)
			subprocess.call("rm %s" % (os.path.join(outDir_5_indiv_lib, (lib + '_unpaired.fastq.gz'))), shell=True)
			print('Finalized removing reads lib %s' % (lib))
			sys.stdout.flush()
		else:
			pass
	fh_f.close()
	fh_r.close()
	fh_u.close()

# Remove low-complexity reads
outDir_6 = os.path.join(results_path, "6.complex_checked")
if (complete == 1) or (remove_complex == 1):
	# Print the date + time, when merging started
	print('Remove low-complexity reads started at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
	sys.stdout.flush()
	if os.path.isdir(outDir_6):
		pass
	else:
		subprocess.call("mkdir '%s'" % (outDir_6), shell=True)
	lib_stats_out_F = os.path.join(outDir_6, 'post_low_complex_F.csv')
	lib_stats_out_R = os.path.join(outDir_6, 'post_low_complex_R.csv')
	lib_stats_out_U = os.path.join(outDir_6, 'post_low_complex_U.csv')
	fh_f = open(lib_stats_out_F, 'wt')
	fh_r = open(lib_stats_out_R, 'wt')
	fh_u = open(lib_stats_out_U, 'wt')
	writer_f = csv.writer(fh_f)
	writer_r = csv.writer(fh_r)
	writer_u = csv.writer(fh_u)
	writer_f.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_r.writerow(('Libary', 'Number of reads', 'Average length per read'))
	writer_u.writerow(('Libary', 'Number of reads', 'Average length per read'))
	indiv_done = []
	for lib in libs:
		indiv = lib.rsplit('_')[0]
		outDir_6_indiv = os.path.join(outDir_6, indiv)
		if indiv not in indiv_done:
			subprocess.call("mkdir '%s'" % (outDir_6_indiv), shell=True)
			indiv_done.append(indiv)
		else:
			pass
		outDir_6_indiv_lib = os.path.join(outDir_6_indiv, lib)
		if os.path.isdir(outDir_6_indiv_lib):
			sys.exit("The lib specific folder %s already exists, please specify a new folder or remove old one first" % (outDir_6_indiv_lib))
		else:
			subprocess.call("mkdir '%s'" % (outDir_6_indiv_lib), shell=True)
		os.chdir(outDir_6_indiv_lib)
		lib_R1_P_path_gz = os.path.join(outDir_5, indiv, lib, (lib + '_1_p.fastq.gz'))
		lib_R2_P_path_gz = os.path.join(outDir_5, indiv, lib, (lib + '_2_p.fastq.gz'))
		lib_U_path_gz = os.path.join(outDir_5, indiv, lib, (lib + '_u.fastq.gz'))
		# Remove low complexity reads for PE and SE
		low_complex_1 = flagLowComplexity(lib_R1_P_path_gz, cut_off_complexity)
		low_complex_2 = flagLowComplexity(lib_R2_P_path_gz, cut_off_complexity)
		low_complex_u = flagLowComplexity(lib_U_path_gz, cut_off_complexity)
		lowcomp_both = low_complex_1 & low_complex_2
		lowcomp_1_only = low_complex_1 - low_complex_2
		lowcomp_2_only = low_complex_2 - low_complex_1
		# New FastQ files
		out1_P = open('%s_R1.fastq' % lib, 'w')
		out2_P = open('%s_R2.fastq' % lib, 'w')
		out_U_new = open('%s_U.fastq' % lib, 'w')
		out_low_complex = open('%s_low_complex.fastq' % lib, 'w')
		# Open Old FastQ files
		in1_P = gzip.open(lib_R1_P_path_gz, 'r')
		in2_P = gzip.open(lib_R2_P_path_gz, 'r')
		inU = gzip.open(lib_U_path_gz, 'r')
		# Now go over the old file and remove low duplicate sequences
		while 1:
			# Store read data
			nameline_1 = in1_P.readline()
			seq_1 = in1_P.readline()
			del_1 = in1_P.readline()
			qual_1 = in1_P.readline()
			nameline_2 = in2_P.readline()
			seq_2 = in2_P.readline()
			del_2 = in2_P.readline()
			qual_2 = in2_P.readline()
			if not nameline_1: break
			rname_1 = nameline_1.strip().split(' ')[0]
			rname_2 = nameline_2.strip().split(' ')[0]
			# Compare to list
			if rname_1 in lowcomp_both:
				#continue
				out_low_complex.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
				out_low_complex.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
			elif rname_1 in lowcomp_1_only:
				out_U_new.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
				out_low_complex.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
			elif rname_2 in lowcomp_2_only:
				out_U_new.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
				out_low_complex.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
			else:
				out1_P.write('%s\n%s%s%s' % (rname_1, seq_1, del_1, qual_1))
				out2_P.write('%s\n%s%s%s' % (rname_2, seq_2, del_2, qual_2))
		while 1:
			nameline_u = inU.readline()
			seq_u = inU.readline()
			del_u = inU.readline()
			qual_u = inU.readline()
			if not nameline_u: break
			rname_u = nameline_u.strip().split(' ')[0]
			if rname_u in low_complex_u:
				out_low_complex.write('%s\n%s%s%s' % (rname_u, seq_u, del_u, qual_u))
			else:
				out_U_new.write('%s\n%s%s%s' % (rname_u, seq_u, del_u, qual_u))
		in1_P.close()
		in2_P.close()
		inU.close()
		out1_P.close()
		out2_P.close()
		out_U_new.close()
		out_low_complex.close()
		gzip_fn(pigz, threads, os.path.join(outDir_6_indiv_lib, (lib + '_R1.fastq')))
		gzip_fn(pigz, threads, os.path.join(outDir_6_indiv_lib, (lib + '_R2.fastq')))
		gzip_fn(pigz, threads, os.path.join(outDir_6_indiv_lib, (lib + '_U.fastq')))
		gzip_fn(pigz, threads, os.path.join(outDir_6_indiv_lib, (lib + '_low_complex.fastq')))
		x, y = getReadStats(os.path.join(outDir_6_indiv_lib, (lib + '_R1.fastq.gz')))
		writer_f.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_6_indiv_lib, (lib + '_R2.fastq.gz')))
		writer_r.writerow((lib, x, y))
		x, y = getReadStats(os.path.join(outDir_6_indiv_lib, (lib + '_U.fastq.gz')))
		writer_u.writerow((lib, x, y))
		if limit_storage == 1:
			print('Start removing reads lib %s' % (lib))
			sys.stdout.flush()
			subprocess.call("rm %s" % (lib_R1_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_R2_P_path_gz), shell=True)
			subprocess.call("rm %s" % (lib_U_path_gz), shell=True)
			print('Finalized removing reads lib %s' % (lib))
			sys.stdout.flush()
		else:
			pass
	fh_f.close()
	fh_r.close()
	fh_u.close()

# Print the date + time, when analysis finalized
print('Clean up finished at: {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
