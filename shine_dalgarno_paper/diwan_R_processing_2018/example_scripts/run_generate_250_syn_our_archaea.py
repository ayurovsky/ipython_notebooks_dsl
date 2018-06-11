import sys
import csv
import subprocess
from subprocess import Popen, PIPE, STDOUT

start = int(sys.argv[1])

first_idx = start * 6 + 1 
end_idx = first_idx + 5
print("will generate files from ", first_idx, " to ", end_idx)


input_csv = "Our_NOT_IN_DIWAN_SD_Archaea_WithAllStats_withTompa_v4.csv"
count = -1 
with open(input_csv, 'r') as in_file:
	csvReader = csv.reader(in_file, delimiter=',')
	for row in csvReader:
		count += 1
		if (row[0] == "name"):	
			continue
		if (count < first_idx):
			continue
		if (count > end_idx):
			break
		in_file = "/home/moamin/ncbi_archaea_db/GbArchaea_FNA/" + row[1] + "_cds_from_genomic.fna.gz"
		print(in_file)
		out_file = row[1] + "_random250_syn_shuff_mean_freq.csv"
		try:
			subprocess.check_output(["Rscript", "generate_250_syn.Rscript", in_file, out_file], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError as exc:
			print("Problem executing generate_250_syn.Rscript for ", out_file, "; ", exc.returncode, exc.output) 
		print("finished " + out_file)
			
