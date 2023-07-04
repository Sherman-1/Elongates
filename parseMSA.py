import  os
import subprocess
import yaml
from tqdm import tqdm
import argparse

from utils import process_multiple_records, read_multifasta, remove_temp_file, write_dicts_to_csv
from utils import Cluster 



def analyseClusters(cov, verbose = False, purge = False):

    # Load yaml file
    with open('env.yaml') as f:
        # use safe_load instead load
        dataMap = yaml.safe_load(f)

    # Load regex dict for current organisms, hard coded for now
    re_dict = dataMap["Regex"]["Scer"] 

    cluster_dir = f"work/{cov}/clusters"

    if purge : os.remove(f"output/{cov}_elongates.csv"), os.remove(f"output/{cov}_events.csv")

    for entry in os.scandir(cluster_dir):

        size = entry.name

        if size != "size_0" and size != "size_1":

            cluster_mem = set()
            print(size)
            size_dir = f"work/{cov}/clusters/{size}"
            for cluster_name in tqdm(os.listdir(size_dir)) if verbose else os.listdir(size_dir): # get every cluster for a given cluster size

                bashCommand = f"muscle -quiet -in work/{cov}/clusters/{size}/{cluster_name} -fasta -out work/{cov}/clusters/{size}/temp.fasta"
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                if error != None:
                    print(f"Muscle error : {error}")
                    break

                if cluster_name not in cluster_mem:
                    cluster_mem.add(cluster_name)
                else : 
                    print(f"{cluster_name} has already been seen ... ")

                records = read_multifasta(f"work/{cov}/clusters/{size}/temp.fasta")
                remove_temp_file(cov,size, skip_check=False)
                elongates, events = process_multiple_records(records, cluster_name, size, re_dict)
                break

            break

        break
                # write_dicts_to_csv(elongates, f"output/{cov}_elongates.csv")
                # write_dicts_to_csv(events, f"output/{cov}_events.csv")

if __name__ == "__main__": 

    parser = argparse.ArgumentParser()
    parser.add_argument('--cov', type=str, help='Coverage')
    parser.add_argument('--verbose', action='store_true', help='Let TQDM display progress bar')    
    parser.add_argument('--purge', action='store_true', help = 'Clean the output directory')
    args = parser.parse_args()
    analyseClusters(args.cov, args.verbose, args.purge)