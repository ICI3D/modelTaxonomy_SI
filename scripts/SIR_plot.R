# SIR_plot.R

plotSIR <- function(ts, plotLab, lwd = 3, col = 1, mtime = MAXTIME, minf = MAXINF, type = 'l'){
  par(lwd = lwd, mar = c(5, 5.1, 1, 1), cex.lab = 3, xaxt = 'n', yaxt = 'n')
  plot(ts$time,               # Time on the x axis
       ts$I,                  # Number infected (I) on the y axis
       xlab = "Time",         # Label the x axis
       ylab = NA,             # Label the y axis
       xlim = c(0, mtime),    # Set limits of x axis    
       ylim = c(0, minf),     # Set limits of y axis        
       type = type,            # Use a line plot for continuous time; step plot for discrete time
       bty = "n",             # Remove the box around the plot
       col = col)             # Line color
  axis(1, at = seq(0, mtime, mtime / 3), labels = NA, xaxt = 's', lwd = 2)
  axis(2, at = seq(0, minf, minf / 10), labels = NA, yaxt = 's', lwd = 2)
  mtext(text = plotLab, side = 3, outer = F, adj = -.35, cex = 3, padj = 0.6)
  mtext(text = 'Infected', side = 2, cex = 3, line = 2) # Label the y axis
}
