# GC Content and Sequence Length Analyzer

# Load required package
library(stringr)

# Read FASTA file
file_name <- readline(prompt="Enter FASTA file name: ")

if (!file.exists(file_name)) {
  stop("File does not exist. Please check the file name.")
}

lines <- readLines(file_name)

# Extract headers (lines starting with >)
headers <- lines[grepl("^>", lines)]

# Extract sequences (non-header lines)
sequences <- lines[!grepl("^>", lines)]

# Function to calculate GC content
calculate_gc <- function(seq) {
  g_count <- str_count(seq, "G")
  c_count <- str_count(seq, "C")
  total_length <- nchar(seq)
  gc_content <- (g_count + c_count) / total_length
  return(gc_content)
}
calculate_at <- function(seq) {
  a_count <- str_count(seq, "A")
  t_count <- str_count(seq, "T")
  total_length <- nchar(seq)
  at_content <- (a_count + t_count) / total_length
  return(at_content)
}
# Apply function to all sequences
gc_values <- sapply(sequences, calculate_gc)
at_values <- sapply(sequences, calculate_at)

# Create results data frame
results <- data.frame(
  Sequence = sub("^>", "", headers),
  Length = nchar(sequences),
  GC_Content = gc_values,
  AT_Content = at_values
)

# Plot GC content distribution
hist(results$GC_Content,
     main = "GC Content Distribution",
     xlab = "GC Content",
     col = "lightblue")
hist(results$AT_Content,
     main = "AT Content Distribution",
     xlab = "AT Content",
     col = "pink")
