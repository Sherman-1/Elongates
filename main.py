import yaml
import subprocess
from test import test


def main():


    dataMap = yaml.safe_load(open('env.yaml'))
    for cov in dataMap['Coverages']:

        arg = cov

        subprocess.call(["bash", "./clustering.sh", arg])

    

if __name__ == "__main__":
    main()
    test()