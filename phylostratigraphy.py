#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  5 13:32:18 2022

@author: paul.roginski

As input : 
- A two columns csv file with :
    - col1 : query
    - col2 : matching taxids
- A csv file matching every taxid with its lineage
- The taxid of the organism of interest (optionnal)

As output : a csv file with
 - col1 : query
 - col2 : most recent common node across all matching lineages
 - col3 : diverging lineages, children of the most recent common node
 

The script needs about 5GB of RAM to run smoothly.
"""


#import multiprocessing
import argparse
import subprocess
import urllib.request
import os
import csv
import pandas as pd
import sys




if __name__ == '__main__':
    
#    print("phylostratigraphy.py")

    # Arguments parsing
    parser = argparse.ArgumentParser()
    parser.add_argument("-I", "--matchs", help="A two columns csv file with : -col1 : query, -col2 : matching taxids. All matchs for a query are assumed to be consecutive (BLAST-type)")
    parser.add_argument("-L", "--lineages", help="All lineages file (ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz > fullnamelineage.dmp)")
#    parser.add_argument("-T", "--taxid", help="taxid of the species of interest", default= "FALSE")
#    parser.add_argument("-N", "--ncpus", type=int, help="Nunmber of cpus to use", default= 1)
    parser.add_argument("-O", "--out", help="output file", default= "phylostrat.tsv")
    args = parser.parse_args()
    
    # The current location and name of the NCBI lineages file
    sd_lineages_archive = "ftp://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
    sd_lineages_file = "fullnamelineage.dmp"
    
    
    
    
    # First load all lineages
#    lineages_file = '/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/TEST/fullnamelineage.dmp'
    lineages_file = args.lineages
    if not lineages_file :
    
        print("The lineages file {} is not found. Downloading it...".format(sd_lineages_file), end =" ")
        urllib.request.urlretrieve(sd_lineages_archive, 'lineages.tar.gz')
        bashCommand = "tar -xvf lineages.tar.gz"
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()
        lineages_file = sd_lineages_file
        print("Done.")
    
    print("Loading the lineages...", end =" ")
    lineages = pd.read_csv(lineages_file, sep="|", names=["taxid", "name", "lineage"], index_col=False)
    print("Done.")

    print("Formatting the lineages...", end =" ")
    # Turn the lineage column into a column containning list.
    lineages.lineage = lineages.lineage.str.strip("[ \t;]").str.split('; ')
    print("Done.")
    
    print("Dropping missleading lineages...", end =" ")
    miss_lins = ["other entries","Viruses","unclassified entries"]
    mask = lineages.lineage.apply(lambda x: x[0] not in miss_lins)
    lineages = lineages[mask]
    print("Done.")
    
    
    
    
    def mrcn(curr_query, curr_matchs) :
    # Compute the mrcn for the previous query (= current_query)
        
        # Lineages that match the given query
        current_lin = lineages.lineage[ lineages.taxid.isin( list(set(curr_matchs)) ) ].values.tolist()


        # Most recent common node of theses lineages
        try :
            mrcn = os.path.commonprefix(current_lin)[-1] ### SCRIPT's CORE ###
        except IndexError :
            print("IndexError : there is no common node among the lineages for {}.".format(curr_query))
            return(["NA", "NA"])
            
        # Because of lineages ending with "" :
#        if mrcn == "" : 
#            mrcn = os.path.commonprefix(current_lin)[-2]
        
        # Index of this mrcn in the lineages
        mrcn_ind = current_lin[0].index(mrcn)
        
        # List of the uniq children of the mrcn 
        diverging_nodes = list(set([item[ mrcn_ind +1 ] for item in current_lin if len(item) > mrcn_ind +1 ]))
        
        return([mrcn, diverging_nodes])




#    matchs_file = '/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/RHODO/GCA_005059875.1_reduced.out'
#    matchs_file = '/run/user/4766/gvfs/smb-share:server=store.intra.i2bc.paris-saclay.fr,share=equipes/BIM/MEMBERS/paul.roginski/Eukaryotes/PHYLO/SMITTIUM/GCA_001970855.1/GCA_001970855.1_reduced.out'
    matchs_file = args.matchs
    
    # Then, parse the matchs file. 
    # For each query, gather all its matching taxids,
    # Find the mrcn beewten them.
    # Find the uniq children of the mrcn.
    phylo_dic = {}
    current_matchs = []
    current_query = ""
    c = 0
    with open(matchs_file, 'r') as matchs :
        
        for match in matchs :
            
            match = match.strip().split('\t')
        
            # If this is a new query :
            if match[0] != current_query :
               
                # If this is not the first query :
                if current_query != "" :
                    phylo_dic[ current_query ] = mrcn(current_query, current_matchs)
                    
                    c = c+1
                    sys.stdout.write("{} queries treated   \r".format(c))
                    sys.stdout.flush()
                    
                # Reset variable for the new query 
#                print(match[0])                
                current_query = match[0]
                current_matchs = []

            if len(match) == 2 :
                current_matchs.append(int(match[1]))
                
            
#            c = c+1
#            if c > 1000000 : break


        
        
    # Write the output to the output file
    with open(args.out, 'w') as csv_file:
        
        writer = csv.writer(csv_file, delimiter = "\t")
        writer.writerow(['query', 'mrcn', 'diverging nodes'])
    
        for key in sorted(phylo_dic.keys()):
            writer.writerow([key] + phylo_dic[key])
    
    
    

    

    
    
    