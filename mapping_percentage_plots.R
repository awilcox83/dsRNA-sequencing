library(dplyr)
library(readr)
library(ggplot2)

df <- read_delim("mapping_percentages.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

#Pick only the columns we need for analysis
df <- subset(df, select = c("Name", "Total", "AlignedOnce","AlignedMultiple"))

#Calculate the percentage of reads that map to reference genome
df$Percentage <- round((df$AlignedMultiple + df$AlignedOnce) / df$Total * 100, digits=2)

#Plot the percentage of mapping reads for each sample in a bar chart
dfplot = ggplot(data=df, aes(x=Name, y=Percentage, group=1)) + labs(x="Sample",y="Reads that map to reference genome (%)", title="Mapping Percentages") + ylim(0,100)+ geom_bar(stat="identity")

#Save Plot Image as PDF
ggsave("mapping_percentages.pdf")