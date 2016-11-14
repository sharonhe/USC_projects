# UV-induced mutation Count

1. Download Exome Sequencing data (`somatic mutation.maf`) from TCGA at the link `https://gdc-portal.nci.nih.gov/files/4aee3e32-2802-4e1e-8577-d74b414f30f7`.

2. Open `.maf` file with Sublime Text, remove the last tab delimiter in each line.
 ![image](https://cloud.githubusercontent.com/assets/16218822/19753521/2526ab52-9bba-11e6-9ce5-e40caeff101a.png)

3. Find adjacent bases of mutated base (3bp) for the purpose of distinguishing dipyrimidine and nondipyrimidine sites: 
      a) Copy the content of maf file into Excel, make 4 new columns as following and save as .txt file:
      b) Upload this `.txt` file onto galaxy.org. (choose interval as file type, choose matched genome)
      c) Use Fetch Sequences to extract genomic DNA according to the start & end position in this interval file. Download the output data as interval file. Copy the adjacent bases column into the file used for mutation count in the next step.

4. Count the mutation in dipyrimidine and nondipyrimidine sites: Use Matlab script `SKCM_shan` and Matlab input file `SKCM.txt` to generate output file `mutCount.txt`.


