# Load required package.
library(data.table)

args = commandArgs(trailingOnly = TRUE)

# Read in eigenvectors to plot PCA
x <- fread(args[1])

# Remove column 2 (redundant)
x[, V2 := NULL]

# Load in sample key
y <- fread("/scratch.global/haasx092/gbs_nov_2021/211120_snp_calling_results/ItascaC12_Claudia_sample_key.csv")
# Make a column called identity that keeps only the most important part of the sample names (cultivars).
# Note that it truncates FY-C20 to FY because for most other samples it was easiest to use the hyphen to serve as the delineator
#y[, identity := sub("-.*$", "", sample_name)]

# set column names
setnames(x, c("sample_ID", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                "PC11", "PC12", "PC13", "PC14", "PC15", "PC16", "PC17", "PC18", "PC19", "PC20"))

setnames(y, c("sample_ID", "sample_name", "classification"))

# Shorten sample name (change from directory path format) so that it can be merged with data table y
x[, sample_ID := sub("/.+$", "", sample_ID)]

# Read in eigenvalues (to determine % variation explained by each PC)
v <- fread(args[2])

# Calculate percent variation (note: I didn't bother renaming the columns to something informative since there is only one)
percentVar = c(PC1=v[1, V1] / sum(v$V1), PC2=v[2, V1] / sum(v$V1), PC3=v[3, V1] / sum(v$V1), PC4=v[4, V1] / sum(v$V1), PC5=v[5, V1] / sum(v$V1), PC6=v[6, V1] / sum(v$V1), PC7=v[7, V1] / sum(v$V1), PC8=v[8, V1] / sum(v$V1))

# Merge data tables
x[y, on="sample_ID"] -> z

z[, col := "grey"] # set the default color to grey because there are too many to do it manually in a practical way

# Set colors based on the major cultivars in the plot
z[classification == "R", col := "black"]
z[classification == "S", col := "red"]

# Set plotting character for all points
z[, pch := 16]

# Make the plot
plot_pcs <- function(arg1, arg2, arg3, arg4, arg5, arg6, pch, col){
par(mar = c(4, 4, 2, 16))
plot(arg1, arg2, xlab = paste0(arg3, round(percentVar[arg4]*100), "%"),
                     ylab = paste0(arg5, round(percentVar[arg6]*100), "%"),
                     main = "Claudia's samples and Itasca-C12",
                     col = col,
                     pch = pch,
                     cex = 1.5,
                     yaxt = 'n')
axis(2, las = 1)

par(oma = c(0, 0, 0, 0))
legend("topright", inset = c(-0.2,0.3), xpd = TRUE,
                legend=c("Resistant", "Susceptible", "Itasca-C12"),
                pch=c(16, 16, 16),
                col=c("black", "red", "grey"),
                ncol=1,
                cex=1.2)
}

pdf(args[3], height=12, width=16)
z[, plot_pcs(PC1, PC2, "PC1: ", 1, "PC2: ", 2, col = col, pch = pch)]
z[, plot_pcs(PC2, PC3, "PC2: ", 2, "PC3: ", 3, col = col, pch = pch)]
z[, plot_pcs(PC3, PC4, "PC3: ", 3, "PC4: ", 4, col = col, pch = pch)]
z[, plot_pcs(PC4, PC5, "PC4: ", 4, "PC5: ", 5, col = col, pch = pch)]
z[, plot_pcs(PC5, PC6, "PC5: ", 5, "PC6: ", 6, col = col, pch = pch)]
z[, plot_pcs(PC6, PC7, "PC6: ", 6, "PC7: ", 7, col = col, pch = pch)]
z[, plot_pcs(PC7, PC8, "PC7: ", 7, "PC8: ", 8, col = col, pch = pch)]
dev.off()

# Save data
save(v, x, y, z, percentVar, file=args[4])
