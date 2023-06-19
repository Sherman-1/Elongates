#!/usr/bin/python3.8 



import  os
from Bio import SeqIO
import subprocess
import polars as pl
from utils import extract_info, read_multifasta
import yaml


# Load yaml file
with open('/home/simon.herman/Bureau/Gits/ClusterGenes/env.yaml') as f:
    # use safe_load instead load
    dataMap = yaml.safe_load(f)


os.chdir("/home/simon.herman/Bureau/SCERs_pep")


# Load regex dict for current organisms, hard coded for now
re_dict = dataMap["Regex"]["Scer"] 

for cov in ("cov_0.5"):

    topList = []
    print(cov)

    for size in next(os.walk(f"{cov}/clusters"))[1]: # get directories of each cluster sizes


        if size != "size_0" and size != "size_1":

            cluster_mem = set()
            print(size)
            

            if os.path.exists(f"{cov}/clusters/{size}/temp.fasta"):
                rm = f"rm {cov}/clusters/{size}/temp.fasta"
                process = subprocess.Popen(rm.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                if error != None:
                    print(f"rm error : {error}")
                    break

            i = 0
            for file in os.listdir(f"{cov}/clusters/{size}"): # get every cluster for a given cluster size

                bashCommand = f"muscle -quiet -in {cov}/clusters/{size}/{file} -fasta -out {cov}/clusters/{size}/temp.fasta"
                process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                if error != None:
                    print(f"Muscle error : {error}")
                    break
                   
                
                if (i % 1000 == 0):

                    print(f"File number {i} processed ... ")
                i+=1

                if file not in cluster_mem:
                    cluster_mem.add(file)
                else : 
                    print(f"{file} has already been seen ... ")
                
                try:
                    records = read_multifasta(f"/home/simon.herman/Bureau/SCERs_pep/{cov}/clusters/{size}/temp.fasta")

                except FileNotFoundError:
                    
                    print(f"File error : {cov}, {size}, {file}")
                lists = extract_info(records,file, size, re_dict)
                

                for index in range(len(lists)):
                    topList.append(lists[index])

                rm = f"rm {cov}/clusters/{size}/temp.fasta"
                process = subprocess.Popen(rm.split(), stdout=subprocess.PIPE)
                output, error = process.communicate()
                if error != None:
                    print(f"rm error : {error}")
                    break
            
                
    df = pl.DataFrame(topList, orient = "row")
    df.columns = ["cluster size", "num_cluster", "sequence id", "specie",
                  "peptide_sequence", "seq length", "Nter dashes","Cter dashes",
                  "max length Nter", 
                  "max length Cter", 
                  "Nter elongate","Nter gaps","Nter openings","Nter elongate AA",
                  "Cter elongate","Cter gaps","Cter openings","Cter elongate AA"]
    df.write_csv(f"/home/simon.herman/Bureau/SCERs_pep/{cov}_seqs.csv", separator = "\t")
                
print(df) 


