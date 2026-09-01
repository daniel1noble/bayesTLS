# Build the TLS-family hex stickers.
#
#   bayesTLS -> man/figures/logo.png (+ logo@2x.png)      [ships with the pkg]
#   freqTLS  -> data-raw/freqTLS-logo.png (+ @2x)         [handoff to the
#                                                          itchyshin/freqTLS repo]
#
# Shared motif (the family resemblance): five 4PL survival curves fanning
# across a thermal palette — survival holds at cool temperatures, crashes early
# when hot — with the dashed 50% line and t50 (LT50) points the classical
# z / CTmax quantities are derived from. The packages differ ONLY in how
# uncertainty is drawn, which is exactly how they differ statistically:
#   bayesTLS: a mist of individual posterior draws behind each curve
#   freqTLS:  smooth 95% confidence ribbons computed from the same spread
#
# Design (client-selected "v7 hybrid", 2026-09-01 design panel): deep
# slate-to-steel-blue gradient interior for depth, thin near-white keyline so
# the hex silhouette stays crisp on both light and dark backgrounds. No border
# colour competes with the fan.
#
# Composition is done manually with a native graphics clipping path (R >= 4.2,
# ragg): the full-bleed art is rasterised, drawn inside a viewport clipped to
# the exact pointy-top hexagon (hexb.in proportions), then the keyline is
# stroked on top. (hexSticker was not used: it does not clip full-bleed art to
# the hexagon. And the art must be rasterised BEFORE clipping — a ggplot grob
# pushes internal viewports with rectangular clip = "on", which would replace
# the hexagonal clip path.)
#
# Run from the package root:  Rscript data-raw/make_hexsticker.R

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(magick)
  library(ragg)
  library(sysfonts)
  library(showtext)
})

## --- Selected design parameters (v7 hybrid) ----------------------------------
GRAD_TOP <- "#081420"   # gradient, top
GRAD_BOT <- "#14395C"   # gradient, bottom (steel blue)
KEYLINE  <- "#F8F9FA"   # thin near-white border
KEY_LWD  <- 4
TITLE_PT <- 34
TITLE_Y  <- 0.25        # wordmark centre, fraction of height
BAND_LO  <- 0.325       # curve band, fractions of height
BAND_HI  <- 0.94
PAPER    <- "#F8F9FA"
THERMAL  <- c("#4CC9F0", "#4895EF", "#F9C74F", "#F3722C", "#EF233C")

title_family <- "sans"
tryCatch({
  sysfonts::font_add_google("Jost", "hexfont")
  title_family <- "hexfont"
}, error = function(e) message("Google font unavailable; using 'sans'."))
showtext::showtext_auto()
showtext::showtext_opts(dpi = 300)

## --- Hexagon geometry: pointy-top, hexb.in proportions (43.9 x 50.8 mm) ------
H <- 1200; m <- 26
R <- (H - 2 * m) / 2
W <- ceiling(sqrt(3) * R + 2 * m)
cx <- W / 2; cy <- H / 2
hex_x <- c(cx, cx + sqrt(3) * R / 2, cx + sqrt(3) * R / 2,
           cx, cx - sqrt(3) * R / 2, cx - sqrt(3) * R / 2)
hex_y <- c(cy + R, cy + R / 2, cy - R / 2, cy - R, cy - R / 2, cy + R / 2)
hex_npc_x <- hex_x / W
hex_npc_y <- hex_y / H

## --- Shared curve data (fixed seed: both stickers show the same fan) ---------
set.seed(2026)
p4pl <- function(x, low, up, k, mid) low + (up - low) / (1 + exp(k * (x - mid)))
xg   <- seq(0, 10, by = 0.04)
mids <- seq(7.1, 2.6, length.out = 5)
k0   <- 1.55

draws <- do.call(rbind, lapply(seq_along(mids), function(i) {
  do.call(rbind, lapply(1:55, function(d) {
    data.frame(temp = i, draw = d, x = xg,
               y = p4pl(xg, runif(1, 0, 0.035), runif(1, 0.955, 1),
                        max(0.6, rnorm(1, k0, 0.22)),
                        mids[i] + rnorm(1, 0, 0.38)))
  }))
}))
med <- do.call(rbind, lapply(seq_along(mids), function(i) {
  data.frame(temp = i, x = xg, y = p4pl(xg, 0.01, 0.99, k0, mids[i]))
}))
t50 <- data.frame(temp = seq_along(mids), x = mids, y = 0.5)
# freqTLS ribbon: per-x 95% band over the SAME simulated spread as the mist.
rib <- do.call(rbind, lapply(split(draws, draws$temp), function(dd) {
  qs <- do.call(rbind, lapply(split(dd, dd$x), function(xx)
    data.frame(x = xx$x[1], lo = quantile(xx$y, 0.025),
               hi = quantile(xx$y, 0.975))))
  qs$temp <- dd$temp[1]
  qs
}))

sx <- function(x) x / 10 * W
sy <- function(y) (BAND_LO + y * (BAND_HI - BAND_LO)) * H
for (df in c("draws", "med", "t50")) {
  assign(df, transform(get(df), px = sx(x), py = sy(y)))
}
rib <- transform(rib, px = sx(x), plo = sy(lo), phi = sy(hi))

grad_ras <- matrix(grDevices::colorRampPalette(c(GRAD_TOP, GRAD_BOT))(256),
                   ncol = 1)

## --- One sticker -------------------------------------------------------------
make_sticker <- function(pkg_name, mode, out_base) {
  uncertainty <- if (identical(mode, "freq")) {
    geom_ribbon(data = rib, aes(px, ymin = plo, ymax = phi,
                                fill = factor(temp), group = temp),
                alpha = 0.16, colour = NA)
  } else {
    geom_line(data = draws,
              aes(px, py, group = interaction(temp, draw),
                  colour = factor(temp)),
              linewidth = 0.24, alpha = 0.05)
  }

  art <- ggplot() +
    annotation_raster(grad_ras, 0, W, 0, H, interpolate = TRUE) +
    annotate("segment", x = 0, xend = W, y = sy(0.5), yend = sy(0.5),
             colour = PAPER, linewidth = 0.3, linetype = "22", alpha = 0.30) +
    uncertainty +
    geom_line(data = med, aes(px, py, colour = factor(temp)),
              linewidth = 1.05, alpha = 0.95, lineend = "round") +
    geom_point(data = t50, aes(px, py, fill = factor(temp)), shape = 21,
               colour = "white", stroke = 0.45, size = 2.3, alpha = 0.95) +
    annotate("text", x = cx, y = TITLE_Y * H, label = pkg_name,
             family = title_family, fontface = "bold",
             size = TITLE_PT / ggplot2::.pt, colour = PAPER) +
    scale_colour_manual(values = THERMAL, guide = "none") +
    scale_fill_manual(values = THERMAL, guide = "none") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, W)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, H)) +
    theme_void()

  art_png <- tempfile(fileext = ".png")
  ragg::agg_png(art_png, W, H, res = 300, units = "px", background = GRAD_TOP)
  print(art)
  invisible(dev.off())
  art_ras <- as.raster(image_read(art_png))

  master <- tempfile(fileext = ".png")
  ragg::agg_png(master, W, H, res = 300, units = "px",
                background = "transparent")
  grid.newpage()
  pushViewport(viewport(clip = polygonGrob(hex_npc_x, hex_npc_y)))
  grid.raster(art_ras, 0.5, 0.5, width = unit(1, "npc"),
              height = unit(1, "npc"))
  popViewport()
  grid.polygon(hex_npc_x, hex_npc_y,
               gp = gpar(col = KEYLINE, fill = NA, lwd = KEY_LWD,
                         linejoin = "round"))
  invisible(dev.off())

  img <- image_read(master)
  image_write(img, paste0(out_base, "@2x.png"))
  image_write(image_scale(img, "x600"), paste0(out_base, ".png"))
}

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
make_sticker("bayesTLS", "bayes", "man/figures/logo")
make_sticker("freqTLS",  "freq",  "data-raw/freqTLS-logo")
message("Wrote man/figures/logo{,@2x}.png and data-raw/freqTLS-logo{,@2x}.png ",
        "(font: ", title_family, ")")
