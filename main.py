import yaml
import subprocess
from test import test
from . import parseMSA, prepare_blast, parseBlast


def main():


    dataMap = yaml.safe_load(open('env.yaml'))
    for cov in dataMap['Coverages']:

        subprocess.call(["bash", "./clustering.sh", cov])

        parseMSA.analyseClusters(cov, verbose = True, purge = False)
        prepare_blast.prepare_db(cov)
        subprocess.call(["bash", "./blast.sh", cov])
        parseBlast.parseBlast(cov)

        
    

if __name__ == "__main__":
    main()
    test()