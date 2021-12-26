# README for claudia_analysis

In order to make the file with the full paths to the data (`FASTQ` files), I decided to take an iterative approach to specify which samples to include based on my knowledge of what the sample names are.
```bash
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/5-most*fastq.gz >> claudia_analysis_filenames.txt
```
**Note:** Two of the `FASTQ` files (one forward _R1_ and one reverse _R2_) have an underscore between the "5" and "most" in the sample name. The following line of code accounts for this difference.
```bash
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/5_most*fastq.gz >> claudia_analysis_filenames.txt
```

The next group of samples are either "All R" or "All S".
```bash
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/All*fastq.gz >> claudia_analysis_filenames.txt
```
The remainder of the samples begin with either "NK" or "RS" depending on which farm the samples were collected from.
```bash
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/NK*fastq.gz >> claudia_analysis_filenames.txt
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/RS*fastq.gz >> claudia_analysis_filenames.txt
```
At this point, you can check to make sure all of the files we expect to be there are actually there:
```bash
cat claudia_analysis_filenames.txt | wc -l
```
The result is 84. Since we have 42 samples and each sample has 2 files (one forward read and one reverse read), this matches our expectations.

Since these are Itasca-C12 samples, it might also be helpful to have data from Reneth's GWAS project.
```bash
ls /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/IC12*fastq.gz >> claudia_analysis_filenames.txt
```
If you count the files again, you should get 400. We have 200 samples in total (there are 158 Itasca-C12 samples), so everything is in order.
