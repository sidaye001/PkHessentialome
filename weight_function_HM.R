library(ggplot2)
library(cowplot)

# Generate Ri values
Ri <- seq(-0.5, 0.5, by = 0.02)

# Calculate function values
W <- -8/5 * Ri^2+ 2/5 * Ri +1

##########################New weight function which reduce the effect of TTAA at the 99%##################
###Passing points (0.5,0.8),(0,1),(-0.5,0.4), and drop quickly at the 99% tail and passing(-0.5,0)

highlight_points <- data.frame(
  Ri = c(-0.5, 0, 0.5),
  Function = c(0.4, 1, 0.8),
  label = c("(-0.5, 0.4)", "(0, 1)", "(0.5, 0.8)")
)
# Create a data frame
data <- data.frame(Ri = Ri, Function = W)

#Weight Function: -5/8 * Ri^2 + 2/5 * Ri + 1
# Plot using ggplot2
###Passing points (-0.5,0.8),(0,1),(0.5,0.4)
ggplot(data, aes(x = Ri, y = Function)) +
  geom_line(color = "blue", linewidth = 1) +
  labs(title = "",
       x = "Ri",
       y = "Weight")+
  theme_bw()+ylim(0,1.5)+theme(panel.grid = element_blank())+
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"))+
  geom_point(data = highlight_points, aes(x = Ri, y = Function), color = "red", size = 3) +
  geom_text(data = highlight_points, aes(label = label), 
            vjust = 1, size = 4, fontface = "italic", color = "black") 


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

# Plot the original and modified functions
plot(Ri, W$W_quadratic, type = "l", col = "blue", lwd = 2, xlab = "Ri", ylab = "Weight",
     main = " ", ylim=c(0,1.5), xlim=c(-0.6,0.6))
lines(Ri, W$W_with_drop, col = "red", lwd = 2)
legend("topright", legend = c("Quadratic Function", "Modified Function"), col = c("blue", "red"), lty = 1)

#####Weight function in the figure
plot(Ri, W$W_with_drop, type = "l", col = "black", lwd = 6, xlab = "ri", ylab = "Wi",
     main = " ", ylim=c(0,1.5), xlim=c(-0.6,0.6))

##5X5 inches


