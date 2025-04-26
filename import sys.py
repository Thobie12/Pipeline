import sys
import subprocess
def main(a,b,c,d)
    a=sys.argv[1]
    b=sys.argv[2]
    c=sys.argv[3]
    d=sys.argv[4]

subprocess.run(["scontrol show job" + a], shell=True, check=True, capture_output=True, text=True)
output = result.stdout
if result.stdout.strip== "0":
    print("Job is running")
else
    print("Job is not running, checking sacct for exit code")
    subprocess.run(["scacct --job" + a, "--format=JobID,State", "--noheader"], shell=True, check=True, capture_output=True, text=True)
    print(result.stdout.strip())


result = s*3
result2 = result + 

