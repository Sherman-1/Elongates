import  os
import subprocess
import yaml
from tqdm import tqdm

from utils import process_multiple_records, read_multifasta, remove_temp_file
from utils import Cluster 


def analyseClusters(cov):

    topList = []
    # Load yaml file
    with open('env.yaml') as f:
        # use safe_load instead load
        dataMap = yaml.safe_load(f)

    # Load regex dict for current organisms, hard coded for now
    re_dict = dataMap["Regex"]["Scer"] 
    cluster_dir = f"work/{cov}/clusters"

    print(cluster_dir)
    for entry in os.scandir(cluster_dir):

        size = entry.name

        if size != "size_0" and size != "size_1":

            cluster_mem = set()
            print(size)
            
            for cluster_name in tqdm(os.listdir(f"work/{cov}/clusters/{size}")): # get every cluster for a given cluster size

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


                cluster = Cluster(cluster_name,
                                    size)

                records = read_multifasta(f"work/{cov}/clusters/{size}/temp.fasta")
                remove_temp_file(cov,size, skip_check=True)
                elongates, events = process_multiple_records(records, cluster_name, size, re_dict)

                print(elongates[0])
                break
            break

        break


if __name__ == "__main__": 

    import sys
    cov = sys.argv[1]
    analyseClusters(cov)