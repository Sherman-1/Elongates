import yaml
import subprocess


def main():


    dataMap = yaml.safe_load(open('env.yaml'))
    for cov in dataMap['Coverages']:

        clustering = f"./cluster.sh {cov}"
        process = subprocess.Popen(clustering.split(), stdout=subprocess.PIPE)
        output, error = process.communicate()

if __name__ == "__main__":
    main()