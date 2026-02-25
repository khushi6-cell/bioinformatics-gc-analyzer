# GC Content and Sequence Length Analyzer

# Load required package
library(stringr)

# Read FASTA file
lines <- readLines(file.path(getwd(), "sample_sequences.fasta"))

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

# Apply function to all sequences
gc_values <- sapply(sequences, calculate_gc)

# Create results data frame
results <- data.frame(
  Sequence = sub("^>", "", headers),
  Length = nchar(sequences),
  GC_Content = gc_values
)

print(results)

# Plot GC content distribution
hist(results$GC_Content,
     main = "GC Content Distribution",
     xlab = "GC Content",
     col = "lightblue")