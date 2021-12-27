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

Next, you should make a file with the list of directories. This `txt` file will come in handy for future steps of the GBS analysis.
```bash
ls Sample*/ -d | tr -d / > 211226_claudia_analysis_sample_directory_list.txt
```
This next step is necessary because we are working with paired-end reads. We are doing it because the file `211222_reneth_gwas_sample_names_and_numbers.csv` contains 2 lines per sample (one for the forward read and one for the reverse read).
```bash
awk 'FNR%2' 211226_claudia_analysis_sample_names_and_numbers.csv > 211226_claudia_analysis_file_list_every_other.csv
```
_Make sure you open the resulting file using_ `vi` _to manually remove the header line_. Once that is done, we can make symbolic links (symlinks) to point to the data rather than take up disk space by needlessly duplicating the original files. **Note** that when I analyzed the original dataset, I set `n` to start at 73 and then increased that value by 1 with each iteration of the `while` loop. Since this iteration of the analysis only contains the GWAS samples, there are gaps in the sequence of sample numbers, necessitating a different approach. The approach I used involves extracting the sample number (`Snumber`) from the file name and using that rather than relying on counting iterations through the loop.
```bash
# Make symlinks to GBS data
cat 211226_claudia_analysis_file_list_every_other.csv | cut -f 9 -d / \
	| while read i; do
	STEM=$(echo $i | rev | cut -f 3,4,5,6,7,8 -d "_" | rev)
	Snumber=$(echo $i | rev | cut -f 3 -d "_"| rev | sed 's/^S//g')
	n=$(printf "%04d\n" $Snumber)
	echo $STEM
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R1_001.fastq.gz Sample_$n/Sample_${n}_R1.fq.gz
	ln -s /home/jkimball/data_delivery/umgc/2021-q4/211108_A00223_0697_BHNY3NDSX2/Kimball_Project_008/${STEM}_R2_001.fastq.gz Sample_$n/Sample_${n}_R2.fq.gz
	done
```
In the next step, we will move back to the `R` statistical environment to create a sample key.
```R
# Move back to R
library(data.table)

# Read in data
x <- fread("211226_claudia_analysis_sample_names_and_numbers.csv")

# Change column names
setnames(x, c("sample_number", "sample_name"))

# Add leading zeros
x[, sample_number := sprintf("%04d", as.numeric(sample_number))]
# Add "Sample_" to each sample number
x[, sample_number := paste0("Sample_", sample_number)]

# Remove beginning the beginning part of the filename to remove the part of the path that is no longer necessary to keep
x[, sample_name := sub("^.+Project_008/", "", sample_name)]

# Remove trailing part of filenames (sample names)---ultimately, we only need one line per sample, not two (a consequence of having 2 files per sample for paired-end reads)
x[, sample_name := sub("_[R1].+$", "", sample_name)]
x[, sample_name := sub("_[R2].+$", "", sample_name)]

# Retain unique values only
x <- unique(x)

# Save to CSV
write.csv(x, file="211226_claudia_analysis_sample_key.csv", row.names = FALSE, sep=",", quote=FALSE)
```

This next bit is just to check if the symlinks are all correct. You can do it manually (by going into each directory and typing `ls -lh`, paying special attention to the first and last to make sure the sample numbers match the "S number" in the filenames provided by UMGC. This way of checking it just quickly puts the same information into a single file so you can view them all at once. **Note:** I only did the forward reads (_R1_) because if they are correct, the reverse reads (_R2_) will also be correct.
```bash
for i in $(cat 211226_claudia_analysis_sample_directory_list.txt);
do
ls -lh ${i}/${i}_R1.fq.gz >> check_symlinks_full_paths.txt
done
```
