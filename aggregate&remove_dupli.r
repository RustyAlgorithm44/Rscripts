# Check for missing values
missing_values <- any(is.na(m1$Gene.Symbol))

# Check for empty strings
empty_values <- any(!nzchar(m1$Gene.Symbol))

# Print the results
cat("Any missing values in Symbol column:", missing_values, "\n")
cat("Any empty values in Symbol column:", empty_values, "\n")

m1 <- subset(m1, nzchar(Gene.Symbol))

# Convert the Gene.Symbol column to a factor to ensure correct grouping
m1$Gene.Symbol <- as.factor(m1$Gene.Symbol)

# Select only the numeric columns (excluding "ID" and "Gene.Symbol")
numeric_columns <- sapply(m1, is.numeric) & !colnames(m1) %in% c("ID", "Gene.Symbol")

# Use the aggregate function to calculate the mean for each Gene.Symbol group
aggregated_data <- aggregate(m1[, numeric_columns], by = list(m1$Gene.Symbol), mean)

# Convert back to a data frame if needed
aggregated_data <- as.data.frame(aggregated_data)	

# Print the result
View(aggregated_data)