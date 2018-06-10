from __future__ import division
import subprocess
from subprocess import Popen, PIPE, STDOUT
import csv
import scipy.stats
from os.path import join
import glob
import gzip
import re

csv_with_prefixes = 'all_prefixes.csv'
csv_with_AGGAG_stats_in_orf = 'aggag_orf_stats.csv'
bac_input_path = "/home/moamin/ncbi_bacteria_db/GbBac_FNA/"
arc_input_path = "/home/moamin/ncbi_archaea_db/GbArchaea_FNA/"
output_file_path = "./"

with open(csv_with_prefixes, 'r') as in_file:
	with open(csv_with_AGGAG_stats_in_orf, 'w') as out_file:
		header = ["full_file_prefix", "f1_true", "f1_total", "f1_freq", "f2_true", "f2_total", "f2_freq", "f3_true", "f3_total", "f3_freq", "all_AGGAG", "all_5mers", "AGGAG_freq"]
		csvWriter = csv.writer(out_file, delimiter=',')
		csvWriter.writerow(header)
		csvReader = csv.reader(in_file, delimiter=',')
		for row in csvReader:
			file_prefix = row[0]
			
			missing = 0
			cds_file = ""
			bac_cds_file = bac_input_path + file_prefix + "_cds_from_genomic.fna.gz"
			bac_names = glob.glob(bac_cds_file)
			if (len(bac_names) < 1):
				arc_cds_file = arc_input_path + file_prefix + "_cds_from_genomic.fna.gz"
				arc_names = glob.glob(arc_cds_file)
				if (len(arc_names) < 1):
					print(file_prefix)
					missing = 1
				else:
					cds_file = arc_cds_file
			else:
				cds_file = bac_cds_file

			
			if (not missing):
				sequence_dict = dict()
				print(cds_file)
				with gzip.open(cds_file, 'rt') as fin:
					current_cds = ""
					current_full_region = ""
					for line in fin:
						if (line[0]==">"):
							if (current_cds != ""):
								sequence_dict[current_cds] = current_full_region
							current_cds = line.split(' ')[0][1:]
							current_full_region = ""
						else:
							current_full_region += re.sub(r'([^NACGT])', "N", line.rstrip().upper())
					sequence_dict[current_cds] = current_full_region

				all_f1_true = 0
				all_f1_false = 0
				all_f2_true = 0
				all_f2_false = 0
				all_f3_true = 0
				all_f3_false = 0
				all_total_AGGAG = 0
				all_total_5mers = 0
				for cds in sequence_dict:
					seq = sequence_dict[cds]
					f1_true = 0
					f1_false = 0
					f2_true = 0
					f2_false = 0
					f3_true = 0
					f3_false = 0
					total_AGGAG = 0

					total_5mers =  len(seq) - 4
					all_total_5mers += total_5mers
					total_AGGAG = len(re.findall('(?=AGGAG)',seq))
					all_total_AGGAG += total_AGGAG	
					
					#print(seq)
					#print(total_5mers)
					#print(total_AGGAG)

					# frame 1
					for i in range(12,(len(seq)-3),3):
						if (seq[i:i+3] == "ATG"):
							if (seq[i-12:i-7] == "AGGAG"):
								f1_true += 1
								#print("here")
							else:
								f1_false += 1
					# frame 2
					for i in range(13,(len(seq)-3),3):
						if (seq[i:i+3] == "ATG"):
							if (seq[i-12:i-7] == "AGGAG"):
								f2_true += 1
							else:
								f2_false += 1
					# frame 3
					for i in range(14,(len(seq)-3),3):
						if (seq[i:i+3] == "ATG"):
							if (seq[i-12:i-7] == "AGGAG"):
								f3_true += 1
							else:
								f3_false += 1

					all_f1_true += f1_true
					all_f1_false += f1_false
					all_f2_true += f2_true
					all_f2_false += f2_false
					all_f3_true += f3_true
					all_f3_false += f3_false
				#all_stats = [file_prefix, all_f1_true, all_f1_true + all_f1_false, "{0:.8f}".format(all_f1_true/(all_f1_true + all_f1_false)), all_f2_true, all_f2_true + all_f2_false, "{0:.8f}".format(all_f2_true/(all_f2_true + all_f2_false)), all_f3_true, all_f3_true + all_f3_false, "{0:.8f}".format(all_f3_true/(all_f3_true + all_f3_false)), all_total_AGGAG, all_total_5mers, "{0:.8f}".format(all_total_AGGAG/all_total_5mers)]
				all_stats = [file_prefix, all_f1_true, all_f1_true + all_f1_false, all_f1_true/(all_f1_true + all_f1_false), all_f2_true, all_f2_true + all_f2_false, all_f2_true/(all_f2_true + all_f2_false), all_f3_true, all_f3_true + all_f3_false, all_f3_true/(all_f3_true + all_f3_false), all_total_AGGAG, all_total_5mers, all_total_AGGAG/all_total_5mers]
				print(all_total_AGGAG)
				print(all_total_5mers)
				print(all_f1_true)		
				print(all_f1_false)		
				print(all_f2_true)		
				print(all_f2_false)		
				print(all_f3_true)		
				print(all_f3_false)		
			else:
				all_stats = [file_prefix, "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A", "N/A"] 
			csvWriter.writerow(all_stats)

