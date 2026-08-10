# ---------------------------
# 1. Load Necessary R Packages
# ---------------------------
# Load the required packages
library(corrplot)  # For creating correlation plots

# ---------------------------
# 2. Set Working Directory
# ---------------------------

# Use forward slashes to avoid escape character issues
setwd("Your own Path")

# ---------------------------
# 3. Load and Prepare the Data
# ---------------------------

# Read the correlation data from a CSV file
# Ensure that strings are not converted to factors
corrdata <- read.csv("GenesExpression.csv", stringsAsFactors = FALSE)

# Remove any rows with missing values to ensure clean data
corrdata <- na.omit(corrdata)

# ---------------------------
# 4. Compute Correlation Matrix
# ---------------------------

# Calculate the Spearman correlation matrix
# Spearman method is non-parametric and assesses monotonic relationships
correlation_matrix <- cor(corrdata, method = "spearman")

# ---------------------------
# 5. Define a Function to Plot Correlation Plots
# ---------------------------

# This function generates correlation plots with specified parameters
plot_correlation <- function(cor_matrix, method = "circle", order = NULL, 
                             type = "full", diag = TRUE, 
                             add = FALSE, title = NULL) {
  # Parameters:
  # cor_matrix: The correlation matrix to plot
  # method: Visualization method ("circle", "square", "shade", "ellipse", "number")
  # order: Ordering method for variables (e.g., "original", "AOE", "FPC", "hclust")
  # type: Type of plot ("full", "upper", "lower")
  # diag: Whether to display the diagonal
  # add: If TRUE, adds to an existing plot
  # title: Title of the plot
  
  # Define default title if not provided
  if (is.null(title)) {
    title <- paste("Correlation Matrix (Method:", method, ")")
  }
  
  # Generate the correlation plot
  corrplot(cor_matrix, 
           method = method, 
           order = order, 
           type = type, 
           diag = diag, 
           add = add,
           title = title,
           mar = c(0,0,1,0))  # Adjust margins to accommodate title
}

# ---------------------------
# 6. Generate Various Correlation Plots
# ---------------------------

# 6.1 Basic Correlation Plot
plot_correlation(correlation_matrix, method = "circle", title = "Basic Correlation Plot")

# 6.2 Square Method with FPC Ordering and Upper Triangle
plot_correlation(correlation_matrix, 
                method = "square", 
                order = "FPC", 
                type = "upper", 
                diag = FALSE, 
                title = "Square Method with FPC Ordering (Upper Triangle)")

# 6.3 Shade Method with AOE Ordering and All Upper Triangles
plot_correlation(correlation_matrix, 
                method = "shade", 
                order = "AOE", 
                diag = FALSE, 
                title = "Shade Method with AOE Ordering")

# 6.4 Ellipse Method with AOE Ordering and Upper Triangle
plot_correlation(correlation_matrix, 
                method = "ellipse", 
                order = "AOE", 
                type = "upper", 
                title = "Ellipse Method with AOE Ordering (Upper Triangle)")

# 6.5 Mixed Correlation Plot with Hierarchical Clustering Ordering
# Lower triangle as "shade" and upper triangle as "pie"
plot_correlation(correlation_matrix, 
                method = "shade", 
                type = "full", 
                title = "Mixed Correlation Plot with Hierarchical Clustering Ordering")
corrplot.mixed(correlation_matrix, 
               lower = "shade", 
               upper = "pie", 
               order = "hclust", 
               title = "Mixed Correlation Plot (Shade & Pie)")
