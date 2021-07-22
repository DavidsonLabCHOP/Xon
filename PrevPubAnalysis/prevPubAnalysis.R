library(readr)

# The purpose of this script is to extract pValues for splice junctions of interest from \
# a previously published (palacino et al) analysis of the splicing impact of LMI070. 
# The differential expression results used here (provided as FullDifExpResults.RDS) were generated using
# DESeq2 differential expression analysis of splice junctions found in palacino et al., LMI070 treated 
# datasets vs control datasets.

#Get GRCh38 position information for splice events that palacino defined as significant in their supplementary table 2)
GRCh38_positions <- read_delim("./GRCh38_positions_V2.txt", 
                               "-", escape_double = FALSE, col_names = FALSE, 
                               trim_ws = TRUE)

gene <- read_csv("~/Documents/Davidson_Lab/Alex/PrevPubAnalysis/gene.txt", 
                 col_names = FALSE)
chr <- GRCh38_positions$X1
start <- GRCh38_positions$X2
end <- GRCh38_positions$X3

PrevPubHits <- data.frame(gene, chr, start, end)
colnames(PrevPubHits) <- c("gene", "chromosome", "start", "end")
rm(gene, chr, start, end, GRCh38_positions)

# Read in previously generated diff expression hits
DESeq2_DiffExp_output <- readRDS("./FullDifExpResults.rds")

# Split the rownames column to get the chr, start and end positions for each splice junction
end <- sapply(strsplit(as.character(rownames(DESeq2_DiffExp_output)), "-"), "[[", 2)
split1 <- sapply(strsplit(as.character(rownames(DESeq2_DiffExp_output)), "-"), "[[", 1)
chr <- sapply(strsplit(as.character(split1), ":"), "[[", 1)
start <- sapply(strsplit(as.character(split1), ":"), "[[", 2)

# Add chr start and end values to dataframe as columns
DESeq2_DiffExp_output["chr"] <- chr
DESeq2_DiffExp_output["start"] <- start
DESeq2_DiffExp_output["end"] <- end
rm(chr, start, end)

#write a function to perform the searching
find_hit <- function(chr, end, window) {
  chr <- gsub("chr", "", chr)
  end_range_low <- as.numeric(end)-as.numeric(window)
  end_range_high <- as.numeric(end)+as.numeric(window)
  end <- end_range_low:end_range_high
  
  subsetQuery <- DESeq2_DiffExp_output[chr == DESeq2_DiffExp_output$chr,]

  subsetChrEnd <- subsetQuery[subsetQuery$end %in% end,]
  return(subsetChrEnd)
}

# Run the function by supplying the query chromosome, start position and end position
output <- apply(PrevPubHits, 1, function(x) {find_hit(gsub("chr", "", x[2]), x[3], 10)})
output <- as.data.frame(do.call(rbind, output))

write.csv(output, file = "./prev_position_pvalues.csv")
