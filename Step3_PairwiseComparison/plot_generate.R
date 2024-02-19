# Generate the plot for comparison
generate_plot_comparison = function(comOMGs_test, species_x, species_y){
  plot_comparison = ggplot(comOMGs_test, aes(x = factor(X1), y = factor(X2), fill = value)) +
                    geom_tile() + 
                    geom_text(aes(label = value), color = "black") +
                    scale_fill_gradient(low = "white", high = "darkgreen") +
                    geom_rect(data = subset(comOMGs_test, test == "Reject"), 
                    aes(xmin = as.numeric(factor(X1)) - 0.5, xmax = as.numeric(factor(X1)) + 0.5, 
                        ymin = as.numeric(factor(X2)) - 0.5, ymax = as.numeric(factor(X2)) + 0.5), 
                        fill = NA, color = "red", size = 1) +
                    theme(legend.position="right", 
        axis.line.x = element_line(size = 1),
        axis.line.y = element_line(size = 1),
        axis.text.x = element_text(color="black", face = "bold",size=13, angle=90, vjust = 0.4, hjust = 1), 
        axis.text.y = element_text(color="black", face = "bold",size=13, angle=0, vjust = 0.5, hjust = 1),
        axis.ticks.length = unit(0.1,"cm"),
        axis.ticks = element_line(size = 1),
        legend.title = element_text(size=10, face= "bold"),
        legend.text = element_text(size=8),
        plot.title = element_text(size=14, face= "bold", colour= "blue", hjust = 0.5),
        axis.title.x = element_text(size=18, face="bold", colour = "darkgreen", vjust = 1),    
        axis.title.y = element_text(size=18, face="bold", colour = "darkgreen")) +
                    labs(x = species_x, y = species_y)

  return(plot_comparison)
}