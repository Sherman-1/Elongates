import subprocess
import polars as pl 


def test():

    command = 'find . -type f -name "*.pep" -exec awk \'/^>/ {print $0}\' {} \\; | wc -l'
    result = subprocess.run(command, shell=True, capture_output=True, text=True)

    # Check if the command was successful
    if result.returncode != 0:
        # The command failed, the error message is in result.stderr
        print(f"Command failed with error: {result.stderr}")

    data = pl.read_csv("output/yeasts_stats.csv")
    count = 0
    for line in data.iter_rows():
        count += line[0]*line[1]


    if int(count) == int(result.stdout):
        print("Test passed")
        return 0
    else:
        print("Test failed")
        return 1