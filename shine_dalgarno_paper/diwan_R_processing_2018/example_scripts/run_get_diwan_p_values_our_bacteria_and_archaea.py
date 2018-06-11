import sys
import csv
import subprocess
from subprocess import Popen, PIPE, STDOUT

input_csv = "Our_NOT_IN_DIWAN_SD_Bacteria_And_Archaea_WithAllStats_withTompa_v4.csv"
output_csv = "Our_NOT_IN_DIWAN_SD_Bacteria_And_Archaea_WithAllStats_withTompa_v4_with_diwan_p_value.csv"
with open(input_csv, 'r') as in_file:
	with open(output_csv, 'w') as out_file:
			csvReader = csv.reader(in_file, delimiter=',')
			csvWriter = csv.writer(out_file, delimiter=',')
			for row in csvReader:
				if (row[0] == "name"):	
					row.append("pvalshuff (aff5 to highest")
					csvWriter.writerow(row)
					continue
				real_file = row[1] + "_real_mean_freq.csv"
				file_250 = row[1] + "_random250_syn_shuff_mean_freq.csv"
				try:
					ret = subprocess.check_output(["Rscript", "get_p_value_with_india_hexamers.Rscript", "accuccuu_output.txt", "accuccuu_hexamer_anti_sd_affinity.csv", file_250, real_file], stderr=subprocess.STDOUT)
					ret = ret.decode("utf-8")
					row.append(ret.strip().split(" ")[1])
					csvWriter.writerow(row)
					#print("_" + ret.strip().split(" ")[1] + "_")
				except subprocess.CalledProcessError as exc:
					print("Problem executing generate_real.Rscript for ", out_file, "; ", exc.returncode, exc.output) 
			
