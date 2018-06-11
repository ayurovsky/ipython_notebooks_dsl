import sys
import csv
import subprocess
from subprocess import Popen, PIPE, STDOUT

input_csv = "Our_NOT_IN_DIWAN_SD_Bacteria_WithAllStats_withTompa_v4.csv"
with open(input_csv, 'r') as in_file:
	csvReader = csv.reader(in_file, delimiter=',')
	for row in csvReader:
		if (row[0] == "name"):	
			continue
		in_file = "/home/moamin/ncbi_bacteria_db/GbBac_FNA/" + row[1] + "_cds_from_genomic.fna.gz"
		print(in_file)
		out_file = row[1] + "_real_mean_freq.csv"
		try:
			subprocess.check_output(["Rscript", "generate_real.Rscript", in_file, out_file], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError as exc:
			print("Problem executing generate_real.Rscript for ", out_file, "; ", exc.returncode, exc.output) 
		print("finished " + out_file)
			
