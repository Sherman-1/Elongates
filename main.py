import yaml
import subprocess
import argparse


from test import test
from parseMSA import analyseClusters
from parseBlast import parseBlast
from prepare_db import prepare_db


def main(verbose, purge):


    dataMap = yaml.safe_load(open('env.yaml'))
    for cov in dataMap['Coverages']:

        subprocess.call(["bash", "./clustering.sh", cov])

        analyseClusters(cov, verbose, purge)
        prepare_db(cov)
        subprocess.call(["bash", "./blast.sh", cov])
        parseBlast(cov)

        
    

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', action='store_true', help='Let TQDM display progress bar')    
    parser.add_argument('--purge', action='store_true', help = 'Clean the output directory')
    args = parser.parse_args()
    
    main(args.verbose, args.purge)
    test()