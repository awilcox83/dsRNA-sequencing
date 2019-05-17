library(dplyr)
library(readr)
library(ggplot2)

file_names <- NULL
file_names <- list.files(".", "*.tab")

for (i in 1:length(file_names)) {
  
#Import Text File for bam-readcount coverage file 
df <- read_delim(file_names[i], "\t", escape_double = FALSE, trim_ws = TRUE)

#Adjust Column Names
colnames(df) <- c("GB", "pos", "N", "count")

#Pick only the columns we need for analysis
df <- subset(df, select = c("GB", "pos", "count"))

#Generate Plot, using facets to diplay coverage plot for each individual segment (different GB numbers correspond to different genome segments)
df_plot <- ggplot(data=df, aes(x=pos, y=count, group=1)) + labs(x="Position",y="Read depth", title=file_names[i]) + facet_wrap(~ GB, scales = "free_x") + geom_line()

csv_file <- paste(file_names[i], ".csv", sep = "")
pdf_file <- paste(file_names[i], ".pdf", sep = "")

#Write data to CSV
write.csv(df, csv_file, row.names = F)

#Save Plot Image as PDF
ggsave(plot = df_plot, pdf_file, device = pdf)
}