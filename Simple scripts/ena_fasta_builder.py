import requests
import re
import subprocess

DEFAULT_AMOUNT = 50
pattern = re.compile(r'/Assemblies/.+')
url_list = []
path = "https://ftp.ebi.ac.uk/pub/databases/ENA2018-bacteria-661k"
filename = "/sampleid_assembly_paths.txt"
out_filename = "final_file.fa"

with requests.get(path + filename, stream=True, timeout=30) as response:
    response.raise_for_status()

    for i, line in enumerate(response.iter_lines(decode_unicode=True)):

        if i >= DEFAULT_AMOUNT:
            break
        
        match = re.search(pattern, line.strip())
        if match:
            url_list.append(match.group())

with open(out_filename, "wb") as f_out, subprocess.Popen(["bgzip", "-d"], stdin=subprocess.PIPE, stdout=f_out) as process:
    try:
        for url in url_list:
            with requests.get(path + url, stream=True, timeout=30) as response:
                response.raise_for_status()
                
                for chunk in response.iter_content(chunk_size=65536):
                    if chunk:
                        process.stdin.write(chunk)
    finally:
        process.stdin.close()

    process.wait()
    if process.returncode != 0:
        raise RuntimeError(f"bgzip failed with returncode {process.returncode}")
   
while True:
    answer = input("Do you want to compress the resulting FASTA file? [y/N] ")
    answer = answer.strip().lower()
    if answer in ("y", "yes"):
        subprocess.run(["bgzip", "-k", out_filename], check=True)
        break
    elif answer in ("n", "no", ""):
        break
    else:
        print("Invalid answer")