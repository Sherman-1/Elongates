import subprocess
import polars as pl 


__all__ = ["test"]

def test(cov):

    command = 'find . -type f -name "*_corr.pep" -exec awk \'/^>/ {print $0}\' {} \\; | wc -l'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Check if the command was successful
    if result.returncode != 0:
        # The command failed, the error message is in result.stderr
        print(f"Command failed with error: {result.stderr}")

    data = pl.read_csv(f"output/{cov}/yeasts_stats.csv")
    count = 0
    for line in data.iter_rows():
        count += line[0]*line[1]


    if int(count) == int(result.stdout):
        print("Test passed")
        return 0
    else:
        print("Test failed")
        return 1
    

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--cov', type=str, help='Coverage')
    args = parser.parse_args()

    test(args.cov)