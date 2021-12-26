# README for claudia_analysis

## Directory setup
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

The next step is to move into the `R` statistical environment to go from file names to a workable CSV file that will be used in the next step of the directory structure setup.
```R
# Read in data using the data.table package
library(data.table)
fread("claudia_analysis_filenames.txt", header=F) -> x

# Change the column name from V1 to something more informative (filename)
setnames(x, "filename")
# Add a new column called sample_number. It will initially contain the entire filename, but we will work to retain only the sample number
x[, sample_number := filename]
# Strip off first part of filename until sample number begins (S) but do not include it.
x[, sample_number := sub("^.+[S]", "", sample_number)]
# Strip off end of the filename (after the sample number) ... begins with "_R1" or "_R2"
x[, sample_number := sub("_[R1].+$", "", sample_number)]

# Convert sample numbers to numerical and add leading zeros to all samples (to help with sorting).
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]

# Reorder rows in ascending order
x[order(sample_number)] -> x

# Set column order (my personal preference for sample_number to come first)
setcolorder(x, c("sample_number", "filename")) -> x

# Write output to CSV
write.csv(x, file="211226_claudia_analysis_sample_names_and_numbers.csv", row.names=FALSE, col.names=FALSE, sep=",", quote=FALSE)

# Save table as an R object
save(x, file="211226_claudia_analysis_sample_names_and_numbers.Rdata")
```
After that is done, use the `CSV` file using `bash` to create the directory structure.<br>
**Note:** The `echo $i` part is not really necessary. I just included it to watch the progress.
```bash
cat 211226_claudia_analysis_sample_names_and_numbers.csv | cut -f 1 -d , \
	| while read i; do
	d=Sample_$i
	echo $i
	mkdir -p $d
	done
```
Once that is done, you will probably notice that there is a directory called `Sample_sample_number` which is an artefact of the code. I probably could change the code so that the header isn't interpreted as a sample name, but it's also super easy to just remove it after the code finishes. You can easily remove it with a one-liner:
```bash
rm -rf Sample_sample_number
```
