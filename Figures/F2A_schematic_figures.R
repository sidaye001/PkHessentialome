library(ggplot2)
library(cowplot)

##############fit a sigmoid function to drop quickly at the 99% of the tail
##############fit a sigmoid function to drop quickly at the 99% of the tail
##############fit a sigmoid function to drop quickly at the 99% of the tail
# Sigmoid function for the drop at the 99% tail
sigmoid <- function(x, midpoint, slope) {
  y <- 1 / (1 + exp(-slope * (x - midpoint)))
  return(y)
}

# Quadratic function
quadratic_function <- function(x) {
  return(-8/5 * x^2 + 2/5 * x + 1)
}

# Generate Ri values
Ri <- seq(-0.5, 0.5, by = 0.02)

# Calculate function values for the quadratic part
weight_function <- function(Ri){
  W_quadratic <- quadratic_function(Ri)
  # Define sigmoid parameters for the drop at the beginning
  midpoint_drop <- -0.48  # Adjust as needed
  slope_drop <- -1000     # Adjust as needed for a drop at the beginning
  # Apply sigmoid drop at the beginning
  sigmoid_drop <- sigmoid(Ri, midpoint_drop, slope_drop)
  W_with_drop <- W_quadratic * (1 - sigmoid_drop)
  return(list(W_with_drop=W_with_drop, W_quadratic=W_quadratic))
}

W <- weight_function(Ri)

#####Weight function in the figure
plot(Ri, W$W_with_drop, type = "l", col = "black", lwd = 6, xlab = "ri", ylab = "Wi",
     main = " ", ylim=c(0,1.5), xlim=c(-0.6,0.6))

##5X5 inches

######Bimodal distribution#########
######Bimodal distribution#########
######Bimodal distribution#########
set.seed(123)  # for reproducibility
data <- c(rnorm(1000, mean = -2, sd = 1.5), rnorm(1000, mean = 2, sd = 1.5))
data1 <- rnorm(1000, mean = -2, sd = 1.5)
data2 <- rnorm(1000, mean = 2, sd = 1.5)
# Create a histogram for the bimodal distribution
hist(data, breaks = 30, col = "lightblue", main = "Bimodal Distribution", xlab = "Values", ylab = "Frequency")
# Estimate kernel density
density_est <- density(data)


# Plot the kernel density
plot(density_est, main = " ", col = "black", lwd = 6)

