# Important Libraries
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(tidyverse)
library(DT)
library(skimr)
library(plotly)
library(ggplot2)
library(corrplot)
library(scales)
library(survival)
library(survminer)
library(DESeq2)
library(caret)
library(psych)
library(viridis)
library(scales)
library(caret)
library(glmnet)
library(caTools)
library(caret)
library(randomForest)
library(pheatmap)
setwd("/Users/gayatrimishra/Documents/gaurav_project")

# Downloading the clinical data
gdcProject<- getGDCprojects()
TCGAbiolinks::getProjectSummary("TCGA-BRCA")

query <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab"
)

query_output<- getResults(query) 
GDCdownload(query)

clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)


# Accessing the patient data out of 9 files 
pat<- clinical.BCRtab.all$clinical_patient_brca %>% 
  DT::datatable(options = list(scrollX = TRUE, keys = TRUE))
tcga_brca_patient<- pat$x$data
View(tcga_brca_patient)


##################### Exploring the patient tcga brca data ############################

#.  Summary statistics of all the variable
describe(tcga_brca_patient) 
dim(tcga_brca_patient)

unique(tcga_brca_patient$initial_pathologic_dx_year)


nd1 <- tcga_brca_patient %>%
  mutate(her2_and_cent17_cells_count= as.numeric(ifelse(her2_and_cent17_cells_count == "[Not Available]", NA, her2_and_cent17_cells_count)))

# Now you can get a summary of your variable
summary(nd1$her2_and_cent17_cells_count)
describe(nd1$her2_and_cent17_cells_count)

#removing first two rows from the dataset

tcga_brca_patient <- tcga_brca_patient[-(1:2), ]
View(tcga_brca_patient)


#.  dimension of the dataset
dim(tcga_brca_patient)
skim(tcga_brca_patient)
 
table(tcga_brca_patient$gender)

# considering only female observations and excluding the male breast cancer data
tcga_brca_patient<- tcga_brca_patient %>% filter(gender == "FEMALE")



# diagnosis age
tcga_brca_patient$form_completion_date <- as.Date(tcga_brca_patient$form_completion_date)

tcga_brca_patient$birth_days_to <- as.numeric(gsub('-', '', tcga_brca_patient$birth_days_to))
tcga_brca_patient$age_at_diagnosis <- as.numeric(tcga_brca_patient$form_completion_date - tcga_brca_patient$birth_days_to / 365.25)
tcga_brca_patient$current_age<- tcga_brca_patient$birth_days_to/365.25


# ggplot(tcga_brca_patient, aes(x = current_age)) +
#   geom_histogram(aes(y = ..density..), fill = "lightblue", binwidth = 2, color = "black", alpha = 0.7) +
#   geom_density(aes(y = ..density..), color = "black", size = 1,) +
#   geom_vline(aes(xintercept = median(current_age)), linetype = "dashed", color = "blue", size = 1) +
#   labs(title = "Age at First Diagnosis", x = "age", y = "density") +
#   theme_minimal()



# Define the position for your text annotation
text_x <- 77 # X-coordinate for the text
text_y <- 0.030  # Y-coordinate for the text

# Define where you want the arrow to point
line_x <- 70  # X-coordinate on the density line
line_y <- 0.019 # Y-coordinate on the density line

ggplot(tcga_brca_patient, aes(x = current_age)) +
  geom_histogram(aes(y = ..density..), fill = "lightblue", binwidth = 2, color = "black", alpha = 0.7) +
  geom_density(aes(y = ..density..), color = "black", size = 1) +
  geom_vline(aes(xintercept = median(current_age, na.rm = TRUE)), color = "blue", size = 1, linetype = "dashed") +
  labs(title = "Age at First Diagnosis", x = "age", y = "density") +
  annotate("text", x = text_x, y = text_y, label = "Kernel Density Estimate", size = 3, hjust = 0.5, vjust = 0.5) +
  geom_segment(aes(x = text_x, y = text_y, xend = line_x, yend = line_y), 
               arrow = arrow(type = "closed", length = unit(0.2, "inches")), 
               color = "black") +
  theme_minimal()



summary(tcga_brca_patient$current_age)

library(ggplot2)


# boxplot(tcga_brca_patient$current_age, 
#         main="Box Plot of Age",
#         xlab="Age",
#         ylab="Years",
#         col="lightblue",
#         border="black")





  iqr <- IQR(tcga_brca_patient$current_age, na.rm = TRUE)
  upper_quartile <- quantile(tcga_brca_patient$current_age, 0.75, na.rm = TRUE)
  lower_quartile <- quantile(tcga_brca_patient$current_age, 0.25, na.rm = TRUE)
  upper_whisker <- upper_quartile + 1.5 * iqr
  lower_whisker <- lower_quartile - 1.5 * iqr
  
  tcga_brca_patient$outlier <- with(tcga_brca_patient, ifelse(current_age > upper_whisker | current_age < lower_whisker, as.character(current_age), NA))
  ggplot(tcga_brca_patient, aes(x = factor(1), y = current_age)) +
    geom_boxplot(fill = "lightblue", color = "black", width = 0.5) +
    geom_jitter(aes(color = current_age), width = 0.1, size = 1, alpha = 0.2) +
    labs(title = "Box Plot of Age at First Diagnosis",
         x = "",  # Since it's a dummy variable, we leave the x-axis label blank.
         y = "Age (Years)") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    scale_color_gradient(low = "blue", high = "red") +
    stat_boxplot(geom ='errorbar') +
    geom_text(aes(label = outlier), na.rm = TRUE, position=position_dodge(width=0.75), hjust=-0.3, size=2.5, vjust=0, check_overlap = TRUE)

  
  
  
  
  
  
  
  library(ggplot2)
  
  # Assuming 'tcga_brca_patient' is your dataframe
  # Calculate median and mean
  median_age <- median(tcga_brca_patient$current_age, na.rm = TRUE)
  mean_age <- mean(tcga_brca_patient$current_age, na.rm = TRUE)
  
  # Create the boxplot
 agebox<- ggplot(tcga_brca_patient, aes(x = factor(1), y = current_age)) +
    geom_boxplot(fill = "lightblue", color = "black", width = 0.5) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.2, color = "darkgray") +
    labs(title = "Box Plot of Age at First Diagnosis",
         x = "",  # No x-axis label needed
         y = "Age (Years)") +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(size = 12),
          axis.title.y = element_text(size = 12)) +
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_text(x = 1, y = median_age, label = paste("Median:", round(median_age, 1)), vjust = 2, hjust = -1.5, size = 3.8) +
    geom_text(x = 1, y = mean_age, label = paste("Mean:", round(mean_age, 1)), vjust = -2, hjust = -1.5, size = 3.8, color = "blue") +
    coord_cartesian(ylim = c(lower_quartile - 1.5 * iqr, upper_quartile + 1.5 * iqr))  # Adjust y-axis limits based on IQR
  
 print(agebox)
  ggsave("agebox plot.png")

  
  
  
  
  
  
  
  

#. Table of gender and its counts

table(tcga_brca_patient$gender)

# Pie Chart of Different races
race_counts <- table(tcga_brca_patient$race)
race_data <- data.frame(Race = names(race_counts), Count = as.numeric(race_counts))
pie_chart <- plot_ly(data = race_data, labels = ~Race, values = ~Count, type = "pie")
pie_chart <- pie_chart %>% layout(title = "Pie Chart of Race Counts")
pie_chart


library(ggplot2)
library(dplyr)

# Filter out '[Not Evaluated]' and '[Not Available]' from the race variable
filtered_data <- tcga_brca_patient %>%
  filter(race != "[Not Evaluated]", race != "[Not Available]", race !="AMERICAN INDIAN OR ALASKA NATIVE")

# Create the boxplot with the filtered data
ggplot(filtered_data, aes(x = race, y = current_age)) +
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Box Plot of Current Age by Race",
       x = "Race",
       y = "Current Age (Years)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels if needed

unique(tcga_brca_patient$race)
#.  total counts of cases with tumor or without tumor
result <- tcga_brca_patient %>%
  group_by(tumor_status) %>%
  summarise(total_count = n())
print(result)

#. performing chi-square test between vital status and tumor status
contigency_table<- table(tcga_brca_patient$tumor_status, tcga_brca_patient$vital_status)
chi_square_result<- chisq.test(contigency_table)

print(chi_square_result)

#.  Exploring lymph nodes in the dataset 

# calculating mean, median and range of lymph nodes examined in the diagnosis
summary(tcga_brca_patient$lymph_nodes_examined_count)
tcga_brca_patient$lymph_nodes_examined_count[tcga_brca_patient$lymph_nodes_examined_count == "[Not Available]"] <- NA
tcga_brca_patient$lymph_nodes_examined_count <- as.numeric(tcga_brca_patient$lymph_nodes_examined_count)
mean_value <- mean(tcga_brca_patient$lymph_nodes_examined_count, na.rm = TRUE)
print(mean_value)
median_value <- median(tcga_brca_patient$lymph_nodes_examined_count, na.rm = TRUE)
print(median_value)


range_values <- range(tcga_brca_patient$lymph_nodes_examined_count, na.rm = TRUE)
print(range_values)


table(tcga_brca_patient$lymph_nodes_examined)

#.  lymph nodes examined count in the patients
p <- ggplot(tcga_brca_patient, aes(x = lymph_nodes_examined_count)) +
  geom_histogram(binwidth = 1, fill = "lightblue", color = "black") +
  labs(title = "Lymph Nodes Examined Count")
p <- ggplotly(p)
p

#.  lymph nodes examined count relation with tumor status 
p1 <- ggplot(tcga_brca_patient, aes(x = tumor_status, y = lymph_nodes_examined_count)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Box Plot of Lymph Nodes Examined Count by Tumor Status")
p1<- ggplotly(p1)
p1

#.  lymph nodes examined count relation with vital status
p2 <- ggplot(tcga_brca_patient, aes(x = vital_status, y = lymph_nodes_examined_count)) +
  geom_boxplot(fill = "lightblue", color = "black") +
  labs(title = "Box Plot of Lymph Nodes Examined Count by vital Status")
p2<- ggplotly(p2)
p2
# count of tumors

tumor_pathologic_pt_counts <- table(tcga_brca_patient$ajcc_tumor_pathologic_pt)
tumor_pathologic_pt_counts_df <- as.data.frame(tumor_pathologic_pt_counts)


colnames(tumor_pathologic_pt_counts_df) <- c("AJCC Tumor Pathologic pT", "Count")
print(tumor_pathologic_pt_counts_df)

# tumor type and AJCC cancer staging 
stg<-ggplot(tcga_brca_patient, aes(x = ajcc_staging_edition, fill = ajcc_tumor_pathologic_pt)) +
  geom_bar(position = "fill") +
  labs(title = "Bar Plot of AJCC Staging Edition by Tumor Pathologic pT") +
  theme(legend.title = element_blank())
stg<- ggplotly(stg)
stg


####################. INVASIVE AND NON INVASIVE CATEGORIZATION #############
breast_data <- tcga_brca_patient[c("bcr_patient_barcode", "gender", "race", "ethnicity", "birth_days_to", "vital_status", "death_days_to", "last_contact_days_to", "ajcc_pathologic_tumor_stage", "lymph_nodes_examined", "lymph_nodes_examined_count", "ajcc_tumor_pathologic_pt", "ajcc_nodes_pathologic_pn", "ajcc_metastasis_pathologic_pm","current_age")]
breast_data <- breast_data %>% 
  mutate(vital_status_numeric = ifelse(vital_status == "Alive", 1, 0))
unique(breast_data$vital_status)
View(breast_data)


breast_data$Cancer_Type <- NA  # Initialize the new column
#logic
for (i in 1:nrow(breast_data)) {
  if (!is.na(breast_data$ajcc_tumor_pathologic_pt[i]) &&
      !is.na(breast_data$ajcc_nodes_pathologic_pn[i]) &&
      !is.na(breast_data$ajcc_metastasis_pathologic_pm[i]) &&
      breast_data$ajcc_pathologic_tumor_stage[i] != '[Discrepancy]' &&
      breast_data$ajcc_pathologic_tumor_stage[i] != '[Not Available]') {
    
    if (grepl("Tis|T1a|T1b|T1c", breast_data$ajcc_tumor_pathologic_pt[i]) &&
        grepl("N0 (i-)|N0|i0", breast_data$ajcc_nodes_pathologic_pn[i]) &&
        grepl("M0|M0 (i+)", breast_data$ajcc_metastasis_pathologic_pm[i]) &&
        grepl("Stage 0|Stage I|Stage IA|Stage IB|Stage IIA|Stage IIA1|Stage IIA2|Stage IIB|Stage IIIA|Stage IIIA1|Stage IIIA2|Stage IIIB|Stage IIIC|Stage IIIC1|Stage IIIC2|Stage IIIC3|Stage IV", breast_data$ajcc_pathologic_tumor_stage[i])) {
      breast_data$Cancer_Type[i] <- "Non-Invasive"
    } else {
      breast_data$Cancer_Type[i] <- "Invasive"
    }
  } else {
    breast_data$Cancer_Type[i] <- "Unclassified"  # Handle missing values and specific categories
  }
}
table(breast_data$Cancer_Type)


#visualisation of key variables

#### histogram of Cancer type
bar_plot <- ggplot(breast_data, aes(x = Cancer_Type)) +
  geom_bar(fill = "blue") +
  labs(
    title = "Distribution of Cancer Types",
    x = "Cancer Type",
    y = "Count"
  ) +
  theme_minimal()

# Print the bar plot
print(bar_plot)




type_counts <- table(breast_data$Cancer_Type)


colors <- c("red", "blue")
pie_chart <- plot_ly(labels = names(type_counts), values = type_counts, type = "pie", marker = list(colors = colors)) %>%
  layout(title = "Distribution of Cancer Types")
pie_chart



########### Heatmap of Cancer staging with Tumor ##############
heatmap_data <- table(breast_data$ajcc_pathologic_tumor_stage, breast_data$ajcc_tumor_pathologic_pt)

 heatmap_plot <- ggplot(data = as.data.frame(heatmap_data), aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis(name = "Frequency", option = "C") +
  labs(
    title = "Heatmap of Cancer Staging vs. Pathologic Tumor Staging",
    x = "Pathologic Tumor Stage",
    y = "Cancer Stage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(heatmap_plot)

############ bar graph of pathologic tumor###########

count_data <- as.data.frame(table(breast_data$ajcc_tumor_pathologic_pt))

# # Create a bar chart
# count_data <- count_data %>%
#   mutate(Tumor_Stage_Group = sub("([T][1-4]).*", "\\1", Var1)) %>%
#   group_by(Tumor_Stage_Group) %>%
#   summarize(Count = sum(Freq)) %>%
#   ungroup()
# bar_chart <- ggplot(count_data, aes(x = Tumor_Stage_Group, y = Count)) +
#   geom_bar(stat = "identity", fill = "blue") +
#   labs(
#     title = "Bar Chart of Tumor Stage Counts",
#     x = "Tumor Stage Group",
#     y = "Count"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# print(bar_chart)





breast_data <- breast_data %>%
  mutate(Tumor_Stage_Group = sub("([T][1-4]).*", "\\1", ajcc_tumor_pathologic_pt))

# Now create the boxplot
boxplot <- ggplot(breast_data, aes(x = Tumor_Stage_Group, y = current_age)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Age at Diagnosis by Tumor Stage Group",
    x = "Tumor Stage Group",
    y = "Age at Diagnosis"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(boxplot)


















############### bar graph of pathologic Nodes########
# count_data1 <- as.data.frame(table(breast_data$ajcc_nodes_pathologic_pn))
# count_data1 <- breast_data %>%
#   count(ajcc_nodes_pathologic_pn) %>%
#   mutate(Node_Stage_Group = sub("(N[0-3]).*", "\\1", ajcc_nodes_pathologic_pn),
#          Node_Stage_Group = ifelse(grepl("NX", ajcc_nodes_pathologic_pn), "NX", Node_Stage_Group)) %>%
#   group_by(Node_Stage_Group) %>%
#   summarize(Count = sum(n)) %>%
#   ungroup()
# 
# # Now plot the merged categories
# bar_chart <- ggplot(count_data1, aes(x = Node_Stage_Group, y = Count)) +
#   geom_bar(stat = "identity", fill = "green") +
#   labs(
#     title = "Bar Chart of Node Stage Counts",
#     x = "Node Stage Group",
#     y = "Count"
#   ) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1))
# 
# # Display the bar chart
# print(bar_chart)




breast_data <- breast_data %>%
  mutate(Node_Stage_Group = sub("(N[0-3]).*", "\\1", ajcc_nodes_pathologic_pn),
         Node_Stage_Group = ifelse(grepl("NX", ajcc_nodes_pathologic_pn), "NX", Node_Stage_Group))

# Now create the boxplot
boxplot <- ggplot(breast_data, aes(x = Node_Stage_Group, y = current_age)) +
  geom_boxplot() +
  labs(
    title = "Boxplot of Age at Diagnosis by Node Stage Group",
    x = "Node Stage Group",
    y = "Age at Diagnosis"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the boxplot
print(boxplot)



















################ bar graph of pathologic Metastasis
count_data2 <- as.data.frame(table(breast_data$ajcc_metastasis_pathologic_pm))

# Create a bar chart
bar_chart <- ggplot(count_data2, aes(x = Var1, y = Freq)) +
  geom_bar(stat = "identity", fill = "purple") +
  labs(
    title = "Bar Chart of Metastasis Stage Counts",
    x = "Metastasis Stage",
    y = "Count"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the bar chart
print(bar_chart)






# Create the boxplot without mapping fill to a factor inside aes()
boxplot <- ggplot(breast_data, aes(x = ajcc_metastasis_pathologic_pm, y = current_age)) +
  geom_boxplot(fill = "green") +  # Set the fill color for the boxplot to green directly
  labs(
    title = "Boxplot of Current Age by Metastasis Stage",
    x = "Metastasis Stage",
    y = "Current Age"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) # Rotate x-axis labels

# Display the boxplot
print(boxplot)



















############
View(breast_data)
table(breast_data$ajcc_pathologic_tumor_stage)
breast_data <- breast_data %>%
  mutate(Simplified_Stage = case_when(
    grepl("Stage IA|Stage IB|Stage I$", ajcc_pathologic_tumor_stage) ~ "Stage I",
    grepl("Stage IIA|Stage IIB|Stage II$", ajcc_pathologic_tumor_stage) ~ "Stage II",
    grepl("Stage IIIA|Stage IIIB|Stage IIIC|Stage III$", ajcc_pathologic_tumor_stage) ~ "Stage III",
    grepl("Stage IV", ajcc_pathologic_tumor_stage) ~ "Stage IV",
    grepl("^X$", ajcc_pathologic_tumor_stage) ~ "Stage X",
    TRUE ~ "Unclassified"
  ))
breast_summary <- breast_data %>%
  group_by(Simplified_Stage, Cancer_Type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(Percent = Count / sum(Count))

bar_chart <- ggplot(breast_summary, aes(x = Simplified_Stage, y = Percent, fill = Cancer_Type)) +
  geom_bar(position = "dodge", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Bar Chart of Tumor Stages by Cancer Type",
    x = "Tumo",
    y = "Percentage",
    fill = "Cancer Type"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Display the bar chart
print(bar_chart)




count_table <- breast_data %>%
  count(Simplified_Stage, Cancer_Type) %>%
  group_by(Cancer_Type) %>%
  mutate(Percentage = n / sum(n) * 100) %>%
  ungroup() %>%
  select(Simplified_Stage, Cancer_Type, Count = n, Percentage)

# You can now print this table to the console or write it to a CSV file
print(count_table)












########### Heatmap of Cancer staging with Tumor ##############
breast_data$Tumor_Stage_Group<- count_data$Tumor_Stage_Group
heatmap_data <- table(breast_data$Simplified_Stage, breast_data$Tumor_Stage_Group)

heatmap_plot <- ggplot(data = as.data.frame(heatmap_data), aes(x = Var2, y = Var1, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis(name = "Frequency", option = "C") +
  labs(
    title = "Heatmap of Cancer Staging vs. Pathologic Tumor Staging",
    x = "Pathologic Tumor Stage",
    y = "Cancer Stage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

print(heatmap_plot)




###################
library(dplyr)
library(tidyr)

# Assuming 'data' is your original dataset
# Calculate the count for each combination of cancer stage and pathologic tumor stage
count_data <- breast_data %>%
  count(Simplified_Stage,Tumor_Stage_Group )

# Calculate the total count for normalization
total_count <- sum(count_data$n)

# Add a column with relative frequencies
count_data <- count_data %>%
  mutate(relative_frequency = n / total_count)

# Convert from long to wide format for heatmap plotting
frequency_table <- count_data %>%
  pivot_wider(names_from = Tumor_Stage_Group, values_from = relative_frequency)

# Print the table to the console
print(frequency_table)



cleaned_count_data <- count_data %>%
  group_by(Simplified_Stage) %>%
  mutate(relative_frequency = n / sum(n)) %>%
  ungroup()

# View the cleaned data with relative frequencies
print(cleaned_count_data)


library(dplyr)
library(tidyr)

# Assuming 'count_data' is your dataframe with the structure you've provided
count_data <- count_data %>%
  group_by(Simplified_Stage) %>%
  mutate(relative_frequency = n / sum(n)) %>%
  ungroup()

# If you want to display this in a wide format (which is typical for heatmaps):
frequency_table <- count_data %>%
  pivot_wider(
    names_from = Tumor_Stage_Group,
    values_from = relative_frequency,
    values_fill = list(relative_frequency = 0)  # Fill in NA values with 0
  )

# Display the table
print(frequency_table)

# You can also round the relative frequencies for better readability
frequency_table <- frequency_table %>%
  mutate(across(-Simplified_Stage, round, digits = 2))

print(frequency_table)









total_counts <- count_data %>%
  group_by(Simplified_Stage) %>%
  summarize(total = sum(n))

# Next, we join the total counts back to the original data and calculate relative frequencies
data_with_frequencies <- count_data %>%
  left_join(total_counts, by = "Simplified_Stage") %>%
  mutate(relative_frequency = n / total)

# Now, we pivot the data to get the wide format
frequency_table <- data_with_frequencies %>%
  select(Simplified_Stage, Tumor_Stage_Group, relative_frequency) %>%
  pivot_wider(
    names_from = Tumor_Stage_Group,
    values_from = relative_frequency,
    values_fill = list(relative_frequency = 0) # Fill in NA values with 0
  )


print(frequency_table)




heatmap <- ggplot(melt(frequency_table, id.vars = 'Simplified_Stage'), aes(x = variable, y = Simplified_Stage, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "blue", na.value = "white", name = "Relative\nFrequency") +
  labs(title = "Heatmap of Tumor Stages by Simplified Stage", x = "Pathologic Tumor Stage", y = "Simplified Stage") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the heatmap
print(heatmap)










# If you want to use this table for a heatmap, you can do so with ggplot2
heatmap <- ggplot(count_data, aes(x = Tumor_Stage_Group, y = Simplified_Stage, fill = relative_frequency)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "blue") +
  labs(
    title = "Heatmap of Cancer Staging vs. Pathologic Tumor Staging",
    x = "Pathologic Tumor Stage",
    y = "Cancer Stage",
    fill = "Relative Frequency"
  ) +
  theme_minimal()

# Print the heatmap
print(heatmap)



















#######
table(breast_data$ajcc_nodes_pathologic_pn)
table(breast_data$ajcc_metastasis_pathologic_pm)

########################## SURVIVAL ANALYSIS OF CLINICAL DATA ##############

# Survival Analysis pdf download
clinBRCA <- GDCquery_clinic("TCGA-BRCA", "clinical")
TCGAanalyze_survival(
  data = clinBRCA,
  clusterCol = "gender",
  main = "TCGA Set\n BRCA",
  height = 10,
  width=10
)

# adding a binary column of Cancer Staging
breast_data$Cancer_Type1 <- ifelse(breast_data$Cancer_Type == "Invasive", 1, 0)

####################### Survival analysis based on the Age groups
library(survival)
library(survminer)

# Classifying age into 'young', 'mid-age', and 'old'
breast_data$age_group <- cut(breast_data$current_age, 
                             breaks = c(0, 40, 60, Inf), 
                             labels = c("young", "mid-age", "old"), 
                             include.lowest = TRUE)

breast_data$last_contact_days_to <- gsub("Not Available", NA, breast_data$last_contact_days_to)
breast_data$last_contact_days_to <- as.numeric(breast_data$last_contact_days_to)
breast_data$death_days_to<- as.numeric(breast_data$death_days_to)

# Kaplan-Meier survival analysis stratified by Cancer_Type and age_group
fit <- survfit(Surv(death_days_to, Cancer_Type1) ~  age_group, data = breast_data)

# Plotting the survival curves
g <- ggsurvplot(
  fit, 
  data = breast_data,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal()
)

print(g)

table(breast_data$ajcc_pathologic_tumor_stage)




breast_data$death_years_to <- breast_data$death_days_to / 365.25

# Proceed with your survival analysis and plotting
fit <- survfit(Surv(death_years_to, Cancer_Type1) ~  age_group, data = breast_data)

# Plotting the survival curves in years
g <- ggsurvplot(
  fit, 
  data = breast_data,
  pval = TRUE,
  conf.int = TRUE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time (years)",
  ylab = "Survival probability",
  ggtheme = theme_minimal()
)

print(g)


breast_data1<- breast_data %>% filter(gender == "FEMALE")
breast_data1$tumor_status<- tcga_brca_patient$tumor_status
View(breast_data1)
filtered_breast_data <- breast_data1 %>%
  filter(!(tumor_status == "Tumor Free" & vital_status == "Dead"))
dim(filtered_breast_data)








filtered_breast_data <- breast_data %>%
  filter(ajcc_pathologic_tumor_stage_grouped %in% c("Stage I", "Stage II", "Stage III"))

# Convert time from days to years
filtered_breast_data$death_years_to <- filtered_breast_data$death_days_to / 365.25

# Perform the survival analysis on the filtered data
fit1 <- survfit(Surv(death_years_to, Cancer_Type1) ~ ajcc_pathologic_tumor_stage_grouped, data = filtered_breast_data)

# Plot the survival curves with the time in years
g1 <- ggsurvplot(
  fit1,
  data = filtered_breast_data,
  pval = TRUE,
  conf.int = TRUE,  # Add this argument to show confidence intervals
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time (years)",
  ylab = "Survival probability",
  ggtheme = theme_minimal()
)

# Display the plot with confidence bands
print(g1)












################ Survival Analysis based on the cancer stage 

breast_data <- breast_data %>%
  mutate(ajcc_pathologic_tumor_stage_grouped = case_when(
    ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB") ~ "Stage I",
    ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB") ~ "Stage II",
    ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC") ~ "Stage III",
    ajcc_pathologic_tumor_stage == "Stage IV" ~ "Stage IV",
    TRUE ~ ajcc_pathologic_tumor_stage  # Assumes you want to keep the original value if none of the above conditions are met
  ))
fit1 <- survfit(Surv(death_days_to, Cancer_Type1) ~  ajcc_pathologic_tumor_stage_grouped, data = breast_data)

g1 <- ggsurvplot(
  fit1,
  data = breast_data,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal(),
  
)


print(g1)

############## SURVIVAL ANALYSIS ON OVERALL SURVIVAL ##########
breast_data$overall_survival<- ifelse(breast_data$vital_status == "Alive", 1, 0)

#1
breast_data$death_days_to <- as.numeric(as.character(breast_data$death_days_to))

# Kaplan-Meier survival analysis stratified by Cancer_Type and age_group
fitx<- survfit(Surv(death_days_to, overall_survival) ~  age_group, data = breast_data)

# Plotting the survival curves
a <- ggsurvplot(
  fit, 
  data = breast_data,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal()
)

print(a)

#2
fity <- survfit(Surv(death_days_to, overall_survival) ~  ajcc_pathologic_tumor_stage_grouped, data = breast_data)

b <- ggsurvplot(
  fit1,
  data = breast_data,
  pval = TRUE,
  risk.table = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.y.text = FALSE,
  legend = "right",
  xlab = "Time",
  ylab = "Survival probability",
  ggtheme = theme_minimal(),
  
)


print(b)

######################## LOGISTIC REGRESSION #################

names(breast_data)

data_log <- tcga_brca_patient[c("current_age","ajcc_pathologic_tumor_stage","er_status_ihc_Percent_Positive", "margin_status","tumor_status","vital_status","race","ethnicity","lymph_nodes_examined_count")]
data_log$Cancer_Type1<- breast_data$Cancer_Type1
View(data_log)


model1 <- glm(formula = Cancer_Type1 ~ current_age, family = "binomial", data = data_log)
summary(model1)


model2<- glm(formula = Cancer_Type1 ~ tumor_status, family = "binomial", data = data_log)
summary(model2)

model3<- glm(formula = Cancer_Type1 ~ ajcc_pathologic_tumor_stage, family = "binomial", data = data_log)
summary(model3)

model4<- glm(formula = Cancer_Type1 ~ margin_status, family = "binomial", data = data_log)
summary(model4)

model5<- glm(formula = Cancer_Type1 ~ er_status_ihc_Percent_Positive, family = "binomial", data = data_log)
summary(model5)

model6<- glm(formula = Cancer_Type1 ~ lymph_nodes_examined_count, family = "binomial", data = data_log)
summary(model6)

model7<- glm(formula = Cancer_Type1 ~ race, family = "binomial", data = data_log)
summary(model7)

model8<- glm(formula = Cancer_Type1 ~ ethnicity, family = "binomial", data = data_log)
summary(model8)

model_log<- glm(formula = Cancer_Type1~current_age+
                  tumor_status+
                  ajcc_pathologic_tumor_stage+
                  margin_status+
                  er_status_ihc_Percent_Positive+
                  lymph_nodes_examined_count+
                  race+
                  ethnicity,
                family = "binomial", data= data_log)
summary(model_log)

# logisitic regression model with accuracy 
data_log2 <- breast_data[c("gender", "race","vital_status","ajcc_pathologic_tumor_stage", "lymph_nodes_examined", "lymph_nodes_examined_count", "ajcc_tumor_pathologic_pt", "ajcc_nodes_pathologic_pn", "ajcc_metastasis_pathologic_pm","current_age", "Cancer_Type1")]
str(data_log2)
View(data_log2)
data_log2$lymph_nodes_examined[data_log2$lymph_nodes_examined == "[Not Available]"] <- NA
data_log2$lymph_nodes_examined_count[is.na(data_log2$lymph_nodes_examined_count)] <- 0
data_log2$lymph_nodes_examined <- ifelse(is.na(data_log2$lymph_nodes_examined_count) | data_log2$lymph_nodes_examined_count == 0, "NO", "YES")


ggplot(data_log2, aes(x=lymph_nodes_examined_count,
                      y= Cancer_Type1))+
  geom_jitter(height= 0.5,
              alpha= .1)+
  geom_smooth(method = "glm",
              method.args = list(family= "binomial"))
model<- glm(Cancer_Type1 ~ lymph_nodes_examined_count, data= data_log2 , family="binomial" )
summary(model)

dummy_data <- dummyVars(~., data = data_log2)

encoded_data <- predict(dummy_data, newdata = data_log2)
encoded_data <- as.data.frame(encoded_data)

set.seed(123)  # for reproducibility
split <- sample.split(encoded_data$Cancer_Type1, SplitRatio = 0.7)
train_data <- subset(encoded_data, split == TRUE)
test_data <- subset(encoded_data, split == FALSE)

model_p <- glm(Cancer_Type1 ~ ., data = train_data, family = binomial)
predictions <- predict(model_p, newdata = test_data, type = "response")
predicted_classes <- ifelse(predictions > 0.5, 1, 0)


confusion_matrix <- table(Actual = test_data$Cancer_Type1, Predicted = predicted_classes)
accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
print(confusion_matrix)
print(paste("Accuracy:", accuracy))

View(breast_data)

new_response<-breast_data[c("Cancer_Type1", "bcr_patient_barcode")]
View(new_response)

write.csv(new_response, "new_response.csv", row.names = FALSE)


