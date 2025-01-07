# ShiitakeIGS1

## Introduction

Shiitake mushrooms (Lentinula edodes) have over 200 cultivars and more than 300 IGS1 sequences available in INSD (the International Nucleotide Sequence Databases). This database integrates IGS1 sequences from shiitake strains across China, Japan, and Taiwan. It is designed to assist users in sequencing their shiitake strains and finding the closest matching sequences within the database, facilitating accurate identification and research.

## Installation

### Using the Deployed Version
Access the deployed version directly at [GitHub Pages](https://raingel.github.io/ShiitakeIGS1/pairwise_alignment.html).

### Local Use
To use the database locally:
1. Download this repository to your local computer.
2. Open the `pairwise_alignment.html` file in your web browser.

## Usage

### **✓ Obtaining IGS1 Sequences from Shiitake**

Shiitake DNA can be extracted using standard laboratory methods or commercial kits. Due to the high repeat number in the rDNA region where IGS1 is located, DNA can be successfully amplified from fresh or dried shiitake products, spawn, or spore powder.

![rDNA Diagram](https://raw.githubusercontent.com/Raingel/ShiitakeIGS1/main/assets/rDNA_schematic.jpg)

**1. Recommended Primers and PCR Program for IGS1 Region:**
- **Primers:**
  - LR119: GCGAAGCTATCATCTGCTGGA
  - 5SR2R: TGGTCCCCCACCRTGGTACTAACT
- **PCR Program:**
  1. Initial denaturation: 94°C for 2 minutes
  2. 30-35 cycles of:
     - 94°C for 30 seconds
     - 60°C for 30 seconds
     - 72°C for 60 seconds
  3. Final extension: 72°C for 60 seconds

**2. Recommended Primers and PCR Program for Whole rDNA Region:**
- **Primers:**
  - LR230: CCACAGCCAAGGGAACGGGCTTGG  
*(Note: The published article mistakenly listed this sequence as CCACAGCCAA`GGGGA`CGGGCTTGG.)*  
  - LR220: CCCGCCGTTTACCCGCGCTTGG
- **PCR Program:**
  1. Initial denaturation: 98°C for 2 minutes
  2. 30-35 cycles of:
     - 98°C for 20 seconds
     - 68°C for 10 minutes
  3. Final extension: 68°C for 10 minutes

**Note:**
1. If your shiitake strains have multiple IGS variants within the genome (as most do), you cannot directly send PCR samples for traditional Sanger sequencing. Instead, use a cloning kit to purify or consider using high-throughput sequencing.
2. To amplify the full-length rDNA fragment (approximately 10 kb), it is recommended to use high-efficiency DNA polymerases with proofreading capability.

### **✓ Instructions for Using the Webpage**

1. **Upload Your Sequences**:
   - To upload the IGS1 sequences of your shiitake strains, use the provided input field. Sequences should be formatted in standard FASTA format. Each sequence must begin with a > symbol followed by the sequence title. Multiple sequences can be included by separating them with individual > symbols and titles.

    ![Query Box](https://github.com/Raingel/ShiitakeIGS1/blob/main/assets/query_box.jpg?raw=true)

2. **Align Sequences**:
   - Click on the "Align Sequences" button to perform sequence alignment using the built-in Smith-Waterman algorithm.

3. **View Results**:
   - The results will be displayed in a table showing the most similar cultivars based on the uploaded sequences. The table includes strain names, similarity scores, and other relevant information. You can sort the results by clicking on the column headers.

   - Clicking on the strain name will display detailed alignment results, including the similarity matrix and alignment diagrams.
   - The detailed comparison window provides accession numbers, allowing users to link to the original sequences.

    ![Results Table](https://github.com/Raingel/ShiitakeIGS1/blob/main/assets/result_table.jpg?raw=true)
   - If the user's strain matches a database strain, a typical similarity matrix will look like the one shown below. For example, Query Seq1 may match 22M0001 variant 2, and Query Seq2 may match 22M0001 variant 1.

    ![Detailed Results](https://github.com/Raingel/ShiitakeIGS1/blob/main/assets/detailed.jpg?raw=true)

4. **Demo Mode**:
   - If you want to see a demonstration, click the "demo" button to load default query sequences (two IGS1 variants of MU30990, cultivar 588) and view the process and results.

## Citation
This database is published in *Journal of Basic Microbiology*  
<https://doi.org/10.1002/jobm.202400452>  
If you use these resources, please cite the article.
## Database Updates
The database is designed to accept sequence updates. In the root directory, there is a `shiitake_list.csv` file. Fill in the table with the required information, as shown in the example below, and submit a pull request. After review, the sequences will be added to the database and rebuilt using GitHub Actions.

Example:
| Cultivar                    | Strain                      | Isolation Source | Locality | Vendor                                | Accession 1 | Accession 2 | Accession 3 | Reference               |
|-----------------------------|-----------------------------|------------------|----------|----------------------------------------|-------------|-------------|-------------|-------------------------|
| Akiyama_A221                | Akiyama_A221                | commercial spawn | Japan    | Akiyama Mycological Institute Co. Ltd. | AB251715    | AB251716    |             | (Babasaki et al., 2007) |
| Akiyama_A526                | Akiyama_A526                | commercial spawn | Japan    | Akiyama Mycological Institute Co. Ltd. | AB251721    |             |             | (Babasaki et al., 2007) |
| Akiyama_A567                | Akiyama_A567                | commercial spawn | Japan    | Akiyama Mycological Institute Co. Ltd. | AB251722    |             |             | (Babasaki et al., 2007) |
| Kawamura_shokuyoukin_K5     | Kawamura_shokuyoukin_K5     | commercial spawn | Japan    | Kawamura Shokuyoukin Kenkyujo Co. Ltd. | AB251765    | AB251804    | AB251805    | (Babasaki et al., 2007) |
