# Figure 1 flow diagram (ggplot2, journal style)
library(ggplot2)
td <- tempdir()
png(file.path("W2_results/Fig1_flow_diagram.png"),
    width = 3300, height = 1260, res = 300)
plot.new(); par(mar = c(0, 0, 0, 0))
plot.window(xlim = c(0, 100), ylim = c(0, 40))

b <- function(x, y, w, h, labels, cex = 0.95, border = "#333333", bg = "#ffffff") {
  rect(x, y, x + w, y + h, col = bg, border = border, lwd = 1.3)
  n <- length(labels)
  if (n == 1) {
    text(x + w/2, y + h/2, labels, cex = cex)
  } else {
    for (i in seq_len(n)) text(x + w/2, y + h - (h/(n+1))*i, labels[i], cex = cex)
  }
}
ar <- function(x1, y1, x2, y2) arrows(x1, y1, x2, y2, length = 0.12, lwd = 1.3, angle = 20)
lab <- function(x, y, s, cex = 0.8, col = "#555555") text(x, y, s, cex = cex, col = col, adj = c(0.5, 0))

# Main chain (centered x=30, width 34)
b(13, 33, 34, 6, c("NHANES 2005\u20132010, adults interviewed", "(n = 17,132)"))
ar(30, 33, 30, 29.6); lab(33, 31.3, "aged \u226520 years (n = 17,132)")
b(13, 24.2, 34, 5.4, c("With BSFS bowel phenotype data", "(n = 14,619)"))
ar(30, 24.2, 30, 20.8); lab(33, 22.5, "excluded: missing 25(OH)D (n = 1,238)")
b(13, 15.4, 34, 5.4, c("With serum 25(OH)D", "(n = 13,381)"))
ar(30, 15.4, 30, 12.0); lab(33, 13.7, "excluded: incomplete covariates (n = 115)")
b(13, 6.2, 34, 5.8, c("Analytic sample", "(n = 13,266)"), bg = "#e8f3f1")

# Exclusion side boxes
b(62, 28.6, 32, 4.6, c("Excluded: missing BSFS stool type", "(n = 2,513)"), cex = 0.82, border = "#bbbbbb", bg = "#f7f7f7")
b(62, 19.8, 32, 4.6, c("Excluded: missing serum 25(OH)D", "(n = 1,238)"), cex = 0.82, border = "#bbbbbb", bg = "#f7f7f7")
b(62, 11.0, 32, 4.6, c("Excluded: incomplete covariates", "(n = 115)"), cex = 0.82, border = "#bbbbbb", bg = "#f7f7f7")
segments(47, 27.0, 62, 30.9, lty = 2, col = "#999999", lwd = 1)
segments(47, 18.2, 62, 22.1, lty = 2, col = "#999999", lwd = 1)
segments(47, 9.4, 62, 13.3, lty = 2, col = "#999999", lwd = 1)

# Bottom banner
rect(8, 0.8, 92, 4.2, col = "#eef5f3", border = "#0f766e", lwd = 1.2)
text(50, 2.5, "Constipation-type: 996 (6.80% weighted)   |   Normal: 11,242 (86.56%)   |   Diarrhea-type: 1,028 (6.64%)",
     cex = 0.95, col = "#0b5a54")
dev.off()
cat("Fig1 saved\n")
