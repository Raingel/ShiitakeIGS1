# %%
import pandas as pd
import subprocess
import os
import requests
import xml.etree.ElementTree as ET

def get_gb_by_acc(acc_list, retry=10):
    acc_str = ",".join(acc_list)
    URI = f"https://www.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc_str}&rettype=gb&retmode=xml"

    while retry > 0:
        try:
            print(f"Fetching: {acc_str}")
            r = requests.get(URI)
            root = ET.fromstring(r.text)
            break
        except Exception as e:
            retry -= 1
            print(f"Fetching {acc_str} failed, retry {retry}: {e}")
    else:
        print(f"Fetching {acc_str} failed, skipped")
        return []

    seq_data = []
    for gb_seq in root.findall(".//GBSeq"):
        acc = gb_seq.find(".//GBSeq_accession-version").text
        seq = gb_seq.find(".//GBSeq_sequence").text
        try:
            cultivar = gb_seq.find(".//GBSeq_feature-table//GBFeature_quals//GBQualifier[GBQualifier_name='cultivar']").find("GBQualifier_value").text
        except:
            cultivar = "-"
        seq_data.append({"acc": acc, "seq": seq, "cultivar": cultivar})
    
    return seq_data

def fasta_reader(handle):
    line = handle.readline()
    while line:
        if line[0] == ">":
            title = line[1:].strip()
            seq = ""
            line = handle.readline()
            while line and line[0] != ">":
                seq += line.strip()
                line = handle.readline()
            yield {"title": title, "seq": seq}
        else:
            line = handle.readline()

def _exec(cmd, suppress_output=True):
    print(f"Executing: {cmd}")
    if suppress_output:
        with open(os.devnull, 'w') as DEVNULL:
            out = subprocess.run(cmd, stdout=DEVNULL, stderr=DEVNULL, shell=True)
        return None
    else:
        out = subprocess.run(cmd, stdout=subprocess.PIPE, shell=True)
        print("Output:")
        print(out.stdout.decode('utf-8'))
        print("Exception:")
        print(out.stderr)
        return out.stdout.decode('utf-8'), out.stderr
# %%
ROOT = "./"
raw_df = pd.read_csv(ROOT + "shiitake_list.csv")
raw_df.fillna("-", inplace=True)

os.makedirs(ROOT + "/seqs", exist_ok=True)
seq_pool = open(ROOT + "/seqs/seq_pool.fas", "w")

acc_list = []
seq_info_map = {}
for index, row in raw_df.iterrows():
    acc1 = row["Accession 1"]
    acc2 = row["Accession 2"]
    acc3 = row["Accession 3"]
    seq_title = f'{row["Strain"]}||{row["Cultivar"]}||{row["Locality"]}'
    for i, acc in enumerate([acc1, acc2, acc3], 1):
        if str(acc) != "-":
            acc_list.append(acc)
            seq_info_map[acc] = f"{seq_title}||{i}"

print(f"Total accessions to fetch: {len(acc_list)}")

seq_data_list = get_gb_by_acc(acc_list)

for seq_data in seq_data_list:
    acc = seq_data["acc"]
    #Remove version number
    acc = acc.split(".")[0]
    seq = seq_data["seq"]
    seq_title = seq_info_map[acc]
    seq_pool.write(f'>{seq_title}||{acc}\n{seq}\n')

seq_pool.close()
# %%
mafft_exec = "mafft"
input = os.path.abspath(ROOT + "/seqs/seq_pool.fas")
output = os.path.abspath(ROOT + "/seqs/seq_pool_aln.fas")
mafft_command = f'"{mafft_exec}" --maxiterate 2 "{input}" > "{output}"'
out = subprocess.run(mafft_command, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
print(">>", mafft_command)
print("Output:")
print(out.stdout)
print("Exception:")
print(out.stderr)

msa = open(ROOT + "/seqs/seq_pool_aln.fas", "r")
start, end = 0, 99999
for record in fasta_reader(msa):
    start_tmp = next((i for i, c in enumerate(record["seq"]) if c != "-"), 0)
    end_tmp = next((len(record["seq"]) - 1 - i for i, c in enumerate(reversed(record["seq"])) if c != "-"), len(record["seq"])) + 1
    if start_tmp > start:
        start = start_tmp
    if end_tmp < end:
        end = end_tmp
msa.close()

print(f"Trim the MSA from {start + 1} to {end}")

with open(ROOT + "/seqs/seq_pool_aln.fas", "r") as msa:
    with open(ROOT + "/seqs/seq_pool_aln_trimmed.fas", "w") as msa_trimmed:
        with open(ROOT + "/seqs/seq_pool_aln_trimmed_nogap.fas", "w") as msa_trimmed_nogap:
            for record in fasta_reader(msa):
                msa_trimmed.write(f'>{record["title"]}\n{record["seq"][start:end]}\n')
                msa_trimmed_nogap.write(f'>{record["title"]}\n{record["seq"][start:end].replace("-", "")}\n')

seq_pool = open(ROOT + "/seqs/seq_pool_aln_trimmed_nogap.fas", "r")
seq_dict = {}
for seq in fasta_reader(seq_pool):
    if seq["seq"] in seq_dict:
        seq_dict[seq["seq"]] += f';{seq["title"]}'
    else:
        seq_dict[seq["seq"]] = seq["title"]

with open(ROOT + "/seqs/seq_pool_aln_trimmed_nogap_unique.fas", "w") as handle:
    for key in seq_dict:
        handle.write(f">{seq_dict[key]}\n{key}\n")

#將最新的序列更新到網頁的資料庫內
#序列在網頁內是以javascript的字串儲存，如下
#//db_sequences_start
#        var internalSequences = `>Akiyama_A221||Akiyama_A221||Japan||1||AB251715;Kawamurashiki_O_2||Kawamurashiki_O_2||Japan||1||AB251776;Meiji_369||Meiji_369||Japan||1||AB251751;KRCF829||Akiyama_A800||Japan||1||AB728657;KRCF816||Kinko_241||Japan||2||AB728656;Kinko327||Kinko_327||Japan||1||AB728654;KRCF1101||KRCF1101||Japan||1||AB781593;KRCF603||Mori_290||Japan||1||AB728655;Cr04||Cr04||China||2||MF541556;Cr20||Cr20||China||2||MF541558;Hunong No.1||Hunong No.1||China||1||MF541570;Junxing No.8||Junxing No.8||China||1||MF541574;L04||L04||China||1||MF541578;L241-4||L241-4||China||1||MF541587;L26||L26||China||1||MF541579;L931||L931||China||1||MF541594;Senyuan 8404||Senyuan 8404||China||1||MF541623;Senyuan No.1||Senyuan No.1||China||2||MF541626;Shenxiang No.12||Shenxiang No.12||China||1||MF541634;Shenxiang No.8||Shenxiang No.8||China||1||MF541630;Suxiang No.1||Suxiang No.1||China||1||MF541637;Xiangza No.26||Xiangza No.26||China||2||MF541640;MU30905||320||Pingtung County, Taiwan||1||LC813573
#ttgttctaaagatttgttcaacttgtttgaactttttctttta
#....
#>MU30863||Wild Strain||Nantou County, Taiwan||2||LC813617;MU31006||Wild Strain||-||1||LC813618;MU31003||-||-||2||LC813622;MU31004||-||-||1||LC813623;MU31005||-||-||1||LC813625
#ttgttctaaagatttgttcaaacttgtttgaactttttcttttatttttcttaacggcatgcaccctaagggcaggctttcaatgcattggatgctgtgagttgatattttttgtgactcttttatgatgactctattctcttccccgtcttatacgactatgctgagaaacagaggtaattctgttcgcaacagaatatttgtgtcctttggggatacacttagtgtacttgccagtgccatttagtctgaattctgtcttttgaaaaactactttctgatacactacaactaacaatacttttttctcagaaaccaaacaacaattatggaaagcttctgtctgggtgtgttcaactgtgagttgcattttcttcccaacttttctgcagcttgtgtctatctattcccttatatgatgaaaatttccatacattcgcagtttccacttgagtgtccaattagtactcccctaatgcagaatgttctttatattccccctctataccatgttgaacccataaagcatgtgttgagtgtgggccatctccaaaagagatgcgctggcaggaagcaggggttgacactatgagggttataatgtccttcctttggatgtctgaactactttttcttttcactctcttctttttctcttagagtgctgtaactaattggtcataatcccctccttgcaggtacttatggtatcagtgaagttgtactccggtatatcatcagtaaagtggtgttaatgtgtgttaaattataacagtcagtcagtaaagtgtgttaagttaataacagtccagtcagtaagtgtgttaaattataacagtcagtcagtaagtgtgttaagttaataacagtaagtgtgttgagttaataacagtcagtcagtaagtgtgttaagttaataacagtccagtcagtaaagtgtgttaagtcaataacacagtattataacagtcagtcagtaagtgtgttaagttacacaacagtttggtcagttaagtgtgttaaattatagttactaaagattttgaaaaaagaagggtacatagttgg
#`;
#//db_sequences_end
#將這段字串取代成新的序列
page = "./pairwise_alignment.html"

with open(page, "r") as f:
    lines = f.readlines()
#將原本db_sequences_start到db_sequences_end間的序列取代成新的序列
import re
lines = "".join(lines)
new_seq = open(ROOT + "/seqs/seq_pool_aln_trimmed_nogap_unique.fas", "r").read()
lines = re.sub(r"//db_sequences_start.*//db_sequences_end", f"//db_sequences_start\n        var internalSequences = `{new_seq}`;\n//db_sequences_end", lines, flags=re.DOTALL)
with open(page, "w") as f:
    f.write(lines)
# %%

        

