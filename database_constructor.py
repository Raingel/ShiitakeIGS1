# %%
import pandas as pd
import subprocess
import os
import shutil
import requests
import xml.etree.ElementTree as ET

# %%
def get_gb_by_acc (acc = "LC813555", retry=10):
    #Download the data from genbank
    URI = "https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=" + acc + "&rettype=gb&retmode=xml"

    root = ""
    while retry>0:
        try:
            r = requests.get(URI)
            root = ET.fromstring(r.text)
            break
        except:
            retry = retry - 1
            print(f"Fetching {acc} failed, retry")
    if root == "":
        #Fetch failed
        print(f"Fetching {acc} failed, skipped")
        return None
    #Parse the data
    acc = root.find(".//GBSeq_accession-version").text
    seq = root.find(".//GBSeq_sequence").text
    try:
        cultivar = root.find(".//GBSeq_feature-table//GBFeature_quals//GBQualifier[GBQualifier_name='cultivar']").find("GBQualifier_value").text
    except:
        cultivar = "-"
    return {"acc":acc, "seq":seq, "cultivar":cultivar}

def fasta_reader (handle):
    line = handle.readline()
    while line:
        if line[0] == ">":
            title = line[1:].strip()
            seq = ""
            line = handle.readline()
            while line and line[0] != ">":
                seq += line.strip()
                line = handle.readline()
            yield {"title":title, "seq":seq}
        else:
            line = handle.readline()

def _exec(cmd,suppress_output=True):
    if suppress_output:
        with open(os.devnull, 'w') as DEVNULL:
            out = subprocess.run(cmd, stdout=DEVNULL, stderr=DEVNULL, shell=True)
        return None
    else:
        out = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        print (">>", cmd)
        print("Output:")
        print(out.stdout.decode('utf-8'))
        print("Exception:")
        print(out.stderr)
        return out.stdout.decode('utf-8'), out.stderr

# %%
ROOT = "./"

# %%
raw_df = pd.read_excel(ROOT+"shiitake_list.xlsx")
#Replace NaN to -
raw_df.fillna("-", inplace=True)

# %%
os.makedirs(ROOT+"/seqs", exist_ok=True)
seq_pool = open(ROOT+"/seqs/seq_pool.fas", "w")
for index, row in raw_df[:].iterrows():
    #Get accession No. from col Accession No. 1, Accession No. 2, Accession No. 3
    acc1 = row["Accession 1"]
    acc2 = row["Accession 2"]
    acc3 = row["Accession 3"]
    cultivar = row["Cultivar"]
    seq_title = f'{row["Strain"]}||{row["Cultivar"]}||{row["Locality"]}'
    #Get the data from genbank
    if str(acc1) != "-":
        seq_data = get_gb_by_acc(acc1)
        if seq_data != None:
            seq_pool.write(f'>{seq_title}||1||{acc1}\n{seq_data["seq"]}\n')
    if str(acc2) != "-":
        seq_data = get_gb_by_acc(acc2)
        if seq_data != None:
            seq_pool.write(f'>{seq_title}||2||{acc2}\n{seq_data["seq"]}\n')
    if str(acc3) != "-":
        seq_data = get_gb_by_acc(acc3)
        if seq_data != None:
            seq_pool.write(f'>{seq_title}||3||{acc3}\n{seq_data["seq"]}\n')
seq_pool.close()


# Define the path to the MAFFT executable
#mafft_exec = "./bin/mafft.bat"
#mafft_exec = "./mafft/mafft-win/mafft.bat"
# Turn the exec to absolute path
#mafft_exec = os.path.abspath(mafft_exec)
mafft_exec = "mafft"
input = os.path.abspath(ROOT+"/seqs/seq_pool.fas")
output = os.path.abspath(ROOT+"/seqs/seq_pool_aln.fas")
# Construct the MAFFT command as a single string
mafft_command = f'"{mafft_exec}" --maxiterate 2  "{input}" > "{output}"'
# Run MAFFT
out = subprocess.run(mafft_command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print(">>", mafft_command)
print("Output:")
print(out.stdout)
print("Exception:")
print(out.stderr)

# %%
#Read the MSA and check the shortest sequence, trimmed all the sequences to fit the shortest one
msa = open(ROOT+"/seqs/seq_pool_aln.fas", "r")
start = 0
end = 99999
for record in fasta_reader(msa):
    start_tmp = 0
    end_tmp = -1
    for i in range(len(record["seq"])):
        if record["seq"][i] != "-":
            start_tmp = i
            break
    for i in range(len(record["seq"])):
        if record["seq"][len(record["seq"])-1-i] != "-":
            end_tmp = len(record["seq"])-1-i +1  #Plus 1 is due to python's slice rule
            break
    if start_tmp > start:
        start = start_tmp
    if end_tmp < end:
        end = end_tmp
msa.close()
#Trim the MSA
print(f"Trim the MSA from {start+1} to {end}")
with open(ROOT+"/seqs/seq_pool_aln.fas", "r") as msa:
    with open(ROOT+"/seqs/seq_pool_aln_trimmed.fas", "w") as msa_trimmed:
        with open(ROOT+"/seqs/seq_pool_aln_trimmed_nogap.fas", "w") as msa_trimmed_nogap:
            for record in fasta_reader(msa):
                msa_trimmed.write(f'>{record["title"]}\n{record["seq"][start:end]}\n')
                msa_trimmed_nogap.write(f'>{record["title"]}\n{record["seq"][start:end].replace("-","")}\n')

#read sequnces from seq_pool_aln_trimmed_nogap.fas
seq_pool = open(ROOT+"/seqs/seq_pool_aln_trimmed_nogap.fas", "r")
#Loop over all sequeces, if identical sequences are found, concatenate their titles
seq_dict = {}
for seq in fasta_reader(seq_pool):
    if seq["seq"] in seq_dict:
        seq_dict[seq["seq"]] += f';{seq["title"]}'
    else:
        seq_dict[seq["seq"]] = seq["title"]

with open(ROOT+"/seqs/seq_pool_aln_trimmed_nogap_unique.fas", "w") as handle:
    for key in seq_dict:
        handle.write(f">{seq_dict[key]}\n{key}\n")
# %%
