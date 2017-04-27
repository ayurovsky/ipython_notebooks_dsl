#TODO: to reproduce paper organisms results, only keep ATG in keep_only_valid_orfs

# expected perl: v5
# expected input files:
# *._cds_from_genomic.fna.gz
# keep_only_valid_orfs.pl alisa_efficient_codon_pair_score_hexamers_paper.pl expected_codon_freq.pl
# get_atcg_freq.pl paper_get_5or_6_match_to_sd_pairs.pl paper_get_sd_p_value.pl
# codon_pair.txt codons.txt syn_codons.txt
import subprocess
from subprocess import Popen, PIPE, STDOUT
import csv
import scipy.stats
from os.path import join
import glob

def run_table_generate_pipeline(file_prefix, input_path, output_path, SD): 
    
    # pre-process to get rid of bad orfs
    cds_file = input_path + file_prefix + "_cds_from_genomic.fna.gz"
    processed_cds_file = output_path + file_prefix + "_cds_processed.fa"
    
    # return empty if cds file is missing
    names= glob.glob(cds_file)
    if (len(names) < 1):
        print(cds_file, " is missing, not calculating.")
        return(["-1","-1","-1","-1"])
    try:
        ret = subprocess.check_output(["perl", "keep_only_valid_orfs.pl", cds_file, processed_cds_file], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        print("Problem executing keep_only_valid_orfs for ", cds_file, "; ", exc.returncode, exc.output)
    
    # get the codon pair scores of the cds
    codon_pair_scores_file = output_path + file_prefix + "_pairwise_syn_shuffle_frm123.csv"
    codon_pair_scores_helper_file = output_path + file_prefix + "_obs_exp_frm123.cts"
    try:
        ret = subprocess.check_output(["perl", "alisa_efficient_codon_pair_score_hexamers_paper.pl", "-s", \
                                       processed_cds_file, "-o", file_prefix, "-z", "300"], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as exc:
        print("Problem executing alisa_efficient_codon_pair_score_hexamers for ", cds_file, "; ", exc.returncode, exc.output)
    
    # get the codon bias of the cds
    codon_bias_file = output_path + file_prefix + ".codon_bias"
    try:
        ps = subprocess.Popen("cat " + processed_cds_file + " | perl expected_codon_freq.pl >" + codon_bias_file, \
                              shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        ret = ps.communicate()[0]
        ret = ret.decode("utf-8")
        if (ret != ""):
            print("_", ret, "_")
            print("Problem executing expected_codon_freq.pl for ", cds_file)  
    except:
        print("Problem executing expected_codon_freq.pl for ", cds_file)
    
    # collect the paper stats
    stats = []
    
    # calculate the depletion score
    try:
        ret = subprocess.check_output(["perl", "get_atcg_freq.pl", "-i", processed_cds_file, "-s", SD])
        ret = ret.decode("utf-8")
        try:
            stats.append(str(format(float(ret), '.3f')))
        except:
            print("Problem: get_atcg_freq.pl did not return a float for ", cds_file)
    except subprocess.CalledProcessError as exc:
        print("Problem executing get_atcg_freq.pl for ", cds_file, "; ", exc.returncode, exc.output)
      
    # calculate the codon pair p-value
    patched_SD = SD
    if (len(SD) == 4):
      if (SD[3] == "G"):
          patched_SD = SD + "A"
      else:
          patched_SD = SD = "G"
    codon_pair_scores_only_matches_file = output_path + file_prefix + ".codon_pairs_SD_match"
    try:
        ps = subprocess.Popen("perl paper_get_5or_6_match_to_sd_pairs.pl -s " + patched_SD + " -c " + \
                                codon_pair_scores_file + ">" + codon_pair_scores_only_matches_file, \
                                shell=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        ret = ps.communicate()[0]
        ret = ret.decode("utf-8")
        if (ret != ""):
            print("_", ret, "_")
            print("Problem executing paper_get_5or_6_match_to_sd_pairs.pl for ", cds_file) 
    except:
        print("Problem executing paper_get_5or_6_match_to_sd_pairs.pl for ", cds_file)
    # do the 1 -tailed test p-value calculation (previously done by hand in Excel)
    # read in the array of all codon pair scores minus the ones with stop codons
    all_f1_list = []
    csvfile = open(codon_pair_scores_file)
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        codon1 = row[0][0:3]
        codon2 = row[0][3:6]
        if ((codon1 in ("TAG", "TGA", "TAA")) or (codon2 in ("TAG", "TGA", "TAA"))):
            continue
        all_f1_list.append(float(row[1]))
    #print(all_f1_list)
    #print(len(all_f1_list))
    matches_f1_list = []
    csvfile = open(codon_pair_scores_only_matches_file)
    reader = csv.reader(csvfile, delimiter=',')
    for row in reader:
        matches_f1_list.append(float(row[1]))
    #print(matches_f1_list)
    t_test_stats = scipy.stats.ttest_ind(all_f1_list, matches_f1_list, equal_var = False) # 2-tailed test, get 1 tail by / 2
    codon_pair_p_value = float(t_test_stats.pvalue/2)
    stats.append(str(format(codon_pair_p_value, '.3f')))
    
    # calculate motif depletion p-value
    motif_depletion_p_value = 0.0
    try:
        ret = subprocess.check_output(["perl", "paper_get_sd_p_value.pl", "-s", \
                                       SD, "-f", codon_bias_file], stderr=subprocess.STDOUT)
        motif_depletion_p_value = float(ret.decode("utf-8"))
        try:
            stats.append(str(format(float(motif_depletion_p_value), '.3f')))
        except:
            print("Problem: paper_get_sd_p_value.pl did not return a float for ", cds_file)
    except subprocess.CalledProcessError as exc:
        print("Problem executing alisa_efficient_codon_pair_score_hexamers for ", cds_file, "; ", exc.returncode, exc.output)
    
    # calculate the combined p-value
    obj = scipy.stats.combine_pvalues([codon_pair_p_value, motif_depletion_p_value])
    stats.append(str(format(float(obj[1]), '.3f')))
    
    # output the final file
    #with open(join(output_file_path, file_prefix + ".table3_values"), 'w') as fout:
     #   fout.write(stats)
    
    # cleanup: remove created temporary files
    #TODO optionally: rm processed_cds_file, codon_pair_scores_file, codon_pair_scores_helper_file codon_bias_file \
    # codon_pair_scores_only_matches_file
    # clean up the files
    p = Popen("rm " + processed_cds_file, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    p = Popen("rm " + codon_pair_scores_file, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    p = Popen("rm " + codon_pair_scores_helper_file, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    p = Popen("rm " + codon_bias_file, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    p = Popen("rm " + codon_pair_scores_only_matches_file, shell=True, stdout=PIPE, stdin=PIPE, stderr=STDOUT)
    p.wait()
    
    return (stats)
    

# expected file
csv_with_tompa_sd_added = 'Our_SD_Bacteria_DataFilesInfo_with_Tompa_fixed_up.csv'
csv_with_table_stats_added = 'Our_SD_Bacteria_DataFiles_WithAllStats_with_Tompa.csv'
input_file_path = "/scratch4/moamin/ncbi_bacteria_db/GbBac_FNA/"
output_file_path = "./"


count = 0
with open(csv_with_tompa_sd_added, 'r') as in_file:
    with open(csv_with_table_stats_added, 'w') as out_file: # change to 'a' to append when re-starting
        csvReader = csv.reader(in_file, delimiter=',')
        csvWriter = csv.writer(out_file, delimiter=',')
        for row in csvReader:
            if (row[0] == "name"):
                row = row + ["depletion_score", "codon pair p_value", "codon usage p_value", "combined p_value", "sd_source"]
                csvWriter.writerow(row)
                continue
            count += 1
            #if (count > 2):
            #    break
            #if (count <= 64):
            #    continue
            #if (count > 4):
            #    break
            file_prefix = row[1]
            #if (file_prefix != "GCA_000946735.1_Vibrio_cholerae_G_298_Guinea"):
             #   continue
            #done_SD = "AGGAGG"
            #if ((row[4] != done_SD) and ("-" not in row[4]) and (len(row[4]) ==4)):
            SD = row[4]
            stats = run_table_generate_pipeline(file_prefix, input_file_path, output_file_path, SD)
            row = row + stats + ["from_tompa"]
            print(",".join(row))
            csvWriter.writerow(row)

