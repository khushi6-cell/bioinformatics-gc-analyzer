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

sliding_gc <- function(seq, window_size = 5) {
  
  seq_length <- nchar(seq)
  
  if (window_size > seq_length) {
    stop("Window size is larger than sequence length.")
  }
  
  gc_values <- c()
  
  for (i in 1:(seq_length - window_size + 1)) {
    
    window_seq <- substr(seq, i, i + window_size - 1)
    
    g_count <- str_count(window_seq, "G")
    c_count <- str_count(window_seq, "C")
    
    gc_content <- (g_count + c_count) / window_size
    
    gc_values <- c(gc_values, gc_content)
  }
  
  return(gc_values)
}
# Create results data frame
results <- data.frame(
  Sequence = sub("^>", "", headers),
  Length = nchar(sequences),
  GC_Content = gc_values,
  AT_Content = at_values
)

cat("\nRunning sliding window GC analysis (window size = 5) on first sequence...\n")

first_seq <- sequences[1]

window_gc <- sliding_gc(first_seq, window_size <- as.numeric(readline(prompt="Enter sliding window size: ")))
window_gc <- sliding_gc(first_seq, window_size)

plot(window_gc,
     type = "l",
     col = "blue",
     main = "Sliding Window GC Content",
     xlab = "Window Position",
     ylab = "GC Content")
# Plot GC content distribution
hist(results$GC_Content,
     main = "GC Content Distribution",
     xlab = "GC Content",
     col = "lightblue")
hist(results$AT_Content,
     main = "AT Content Distribution",
     xlab = "AT Content",
     col = "pink")
plot(window_gc,
     type="l",
     col="blue",
     main="Sliding Window GC Content",
     xlab="Window Position",
     ylab="GC Content")
abline(h=mean(window_gc), col="red", lty=2)
