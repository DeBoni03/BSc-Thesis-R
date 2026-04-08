cutoff=function(int,prec){
  d <- data.frame(bin = int,precision = prec)
  l=length(int)
  p1 <- c(d$bin[1], d$precision[1])
  p6 <- c(d$bin[l], d$precision[l])
  m <- (p6[2] - p1[2]) / (p6[1] - p1[1])
  q <- p1[2] - m * p1[1]
  d$diag <- m * d$bin + q
  d$dist <- d$precision - d$diag
  best_cutoff_index <- which.max(d$dist)
  best_cutoff_value <- d$bin[best_cutoff_index]
  print(paste("Suggested CutOff:", best_cutoff_value))
  ggplot(d, aes(x = bin, y = precision)) +
    geom_col(fill = "steelblue", alpha = 0.6) +
    geom_line(color = "blue", size = 1) +
    geom_point(size = 3, color = "blue") +
    geom_line(aes(y = diag), color = "red", linetype = "dashed") +
    geom_segment(aes(xend = bin, yend = diag), color = "gray30") +
    geom_point(data = d[best_cutoff_index, ],aes(y = precision), color = "red", size = 5, shape = 1) +
    annotate("text", x = best_cutoff_value, y = d$precision[best_cutoff_index] + 0.05, label = "Cutoff", color = "red", fontface = "bold") +
    theme_minimal() +
    labs(title="Precision Barplot",y = "Precision", x = "Bin")
}
