# Builds the actiRhythm hex logo (man/figures/logo.png). Run from the package root.

library(ggplot2)
library(showtext)
library(sysfonts)
library(png)

font_add_google("Nunito", "rhythm")
showtext_auto()
showtext_opts(dpi = 500)

blue <- "#236192"; gold <- "#E69F00"; pale <- "#CFE0EC"
size <- 6
baseline <- -0.48

theme_hex <- theme_void() +
  theme(plot.background  = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.margin = margin(0, 0, 0, 0))
coords <- coord_fixed(xlim = c(-1, 1), ylim = c(-1.05, 1.05), expand = FALSE)

# 
ink_extent <- function(label, anchor = -0.9) {
  f <- tempfile(fileext = ".png")
  ggsave(f, ggplot() +
           annotate("text", x = anchor, y = 0, label = label, hjust = 0,
                    size = size, fontface = "bold", family = "rhythm") +
           coords + theme_hex,
         width = 5.08, height = 5.86, units = "cm", dpi = 500, bg = "transparent")
  m <- readPNG(f); w <- dim(m)[2]
  cols <- which(colSums(m[, , 4]) > 0)
  c(left = (min(cols) - w / 2) * 2 / w, right = (max(cols) - w / 2) * 2 / w)
}

# 
full  <- ink_extent("actiRhythm")
acti  <- ink_extent("acti")
actir <- ink_extent("actiR")
shift  <- -(full["left"] + full["right"]) / 2
word_x <- -0.9 + shift
r_x    <- (acti["right"] + actir["right"]) / 2 + shift

hexagon <- data.frame(
  x = c(0, -sqrt(3) / 2, -sqrt(3) / 2, 0, sqrt(3) / 2, sqrt(3) / 2),
  y = c(1, 0.5, -0.5, -1, -0.5, 0.5))
h <- seq(0, 48, length.out = 500)
wave <- data.frame(x = -0.6 + 1.2 * h / 48, y = 0.34 + 0.30 * sin(2 * pi * (h - 6) / 24))
peak <- which.max(wave$y); trough <- which.min(wave$y)

logo <- ggplot() +
  geom_polygon(data = hexagon, aes(x, y), fill = blue, colour = pale, linewidth = 3) +
  annotate("segment", x = -0.66, xend = 0.66, y = 0.34, yend = 0.34,
           colour = pale, linewidth = 0.4, alpha = 0.45) +
  geom_line(data = wave, aes(x, y), colour = "white", linewidth = 2.6, lineend = "round") +
  annotate("point", x = wave$x[peak],   y = wave$y[peak] + 0.07,   colour = gold, size = 4.5) +
  annotate("point", x = wave$x[trough], y = wave$y[trough] - 0.07, colour = pale, size = 3) +
  annotate("text", x = word_x, y = baseline, label = "actiRhythm", colour = "white",
           hjust = 0, size = size, fontface = "bold", family = "rhythm") +
  annotate("text", x = r_x, y = baseline, label = "R", colour = gold,
           hjust = 0.5, size = size, fontface = "bold", family = "rhythm") +
  coords + theme_hex

ggsave("man/figures/logo.png", logo, width = 5.08, height = 5.86,
       units = "cm", dpi = 500, bg = "transparent")
