
data <- read.table("/Users/urileal/ESCI/2nd_year_BDBI/Biophysics/GROUP_PROJECT/ala7.txt", header = TRUE) # Replace "your_text_file.txt" with your file name

write.csv(data, "ala7.csv", row.names = FALSE) 

install.packages("ggplot2") 
library(ggplot2)

data <- read.csv("ala7.csv")

View(data)

# Compute the absolute interaction energy using the provided formula
data$abs_interaction_energy <- with(data, abs(Electric - Electric_Ala) + abs(Vdw - Vdw_Ala) + abs(Solvation - Solvation_Ala) - abs(Solvation_Chain - Solvation_Chain_Ala) - abs(Solvation - Solvation_Ala))

ggplot(data, aes(x = data$Res, y = data$abs_interaction_energy)) +
  geom_bar(stat = "identity", fill = "blue") +
  labs(title = "Analysis of Individual Amino Acid Contributions to Protein Complex", x = "Residues", y = "Change in interaction energy")


threshold_value <- 15  # Adjust this value based on your dataset

# Create a bar plot
ggplot(data, aes(x = Res, y = abs_interaction_energy)) +
  geom_bar(stat = "identity", aes(fill = ifelse(abs_interaction_energy > threshold_value, "High Energy", "Low Energy")), width = 0.5) +
  scale_fill_manual(values = c("High Energy" = "#7FB3D5", "Low Energy" = "blue")) +
  labs(title = "Analysis of Individual Amino Acid Contributions to Protein Complex", x = "Residues", y = "Change in interaction energy") +
  theme_minimal()
