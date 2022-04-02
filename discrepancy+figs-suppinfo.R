library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(schumaker)

# Input sources ----------------------------------------------------------------
handsheet_prctile <- read_csv("cal_flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-prctile-80-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.dat")
handsheet_box <- read_csv("cal_flexnum-prob-0_1-alpha-1-upperfibrethickness-20_biasinvariant-responsevariable-sheetthickness-box-cv-rveerr-0_02-repeats-5_micrometer-FRANK-PTI-analyser-MorFi-fibremorphology-wall-lumen-wall-free_grammage-all_unpressed-unrefined.dat")
exptdata2cal <- read_csv("exptdata2cal_thickness-tissue-ISO-12625-3_grammage-all_unpressed-unrefined.dat")
nobs <- NROW(exptdata2cal)
local_sheet_thickness <- read_csv("local-sheet-thickness-distribution_vs_grammage.dat")

# Select and combine the necessary data for thickness --------------------------
offset <- 1 + nobs * c(
  thickness = 1
)
idx <- 1:nobs

thickness_prctile <- handsheet_prctile[, c(1, offset["thickness"] + idx)]
names(thickness_prctile)[1] <- "grammage"
# Add mean and standard deviation of experimental data at different positions
thickness_prctile <- add_column(
  thickness_prctile,
  experimental = exptdata2cal$thickness_avg,
  .after = "grammage"
)
thickness_prctile <- add_column(
  thickness_prctile,
  intrinsic = exptdata2cal$thickness_correct_avg,
  se = exptdata2cal$thickness_stddev,
  method = "percentile"
)
thickness_box <- handsheet_box[, c(1, offset["thickness"] + idx)]
names(thickness_box)[1] <- "grammage"
thickness_box <- add_column(
  thickness_box,
  experimental = exptdata2cal$thickness_avg,
  .after = "grammage"
)
thickness_box <- add_column(
  thickness_box,
  intrinsic = exptdata2cal$thickness_correct_avg,
  se = exptdata2cal$thickness_stddev,
  method = "box"
)

thickness <- rbind(thickness_prctile, thickness_box, make.row.names = FALSE)

# Leave-one-out cross-validation -----------------------------------------------
# Compute the discrepancy function
SSS_thickness_prctile <- lapply(1:nobs, function(ind) {
  Schumaker(
    thickness_prctile$grammage[-ind],
    t(thickness_prctile$experimental[-ind] - thickness_prctile[-ind, 2 + ind])
  )
})
discrepancy_prctile <- sapply(1:nobs, function(ind) {
  SSS_thickness_prctile[[ind]]$Spline(thickness_prctile$grammage[ind])
})
SSS_thickness_box <- lapply(1:nobs, function(ind) {
  Schumaker(
    thickness_box$grammage[-ind],
    t(thickness_box$experimental[-ind] - thickness_box[-ind, 2 + ind])
  )
})
discrepancy_box <- sapply(1:nobs, function(ind) {
  SSS_thickness_box[[ind]]$Spline(thickness_box$grammage[ind])
})

# Plot settings ----------------------------------------------------------------
fig_width <- 7
fig_asp <- .75
fig_height <- fig_asp * fig_width

# Plot percentiles of local sheet thicknesses vs. grammage ---------------------
# black & white, PDF+EPS
names(local_sheet_thickness) <- c("grammage", "thickness")
df <- dplyr::summarise(
  dplyr::group_by(local_sheet_thickness, grammage),
  prctile = quantile(thickness, probs = c(0.60, 0.70, 0.80, 0.90), type = 8),
  prob = c("60", "70", "80", "90")
)
df$lab <- ""
df$lab[41:44] <- c("60th", "70th", "80th", "90th percentile")

ggplot(df, aes(grammage, prctile, colour = prob, label = lab)) +
  geom_line() +
  geom_text_repel() +
  scale_x_continuous(expand = expansion(add = c(0, 10))) +
  scale_colour_grey(start = 0, end = 0.5) +
  xlab(expression("grammage" ~~ (g / m^2))) +
  ylab(expression(z^s * (phantom() %.% phantom()) ~~ (mu * m))) +
  labs(colour = "percentile") +
  theme(
    legend.position = "none",
    # legend.position = "top",
    # legend.justification = "right"
  )
ggsave("figure4-local-sheet-thickness-percentiles_vs_grammage.pdf",
       width = fig_width, height = 1.25 * fig_height
)
ggsave("figure4-local-sheet-thickness-percentiles_vs_grammage.eps",
       width = fig_width, height = 1.25 * fig_height,
       device = cairo_ps, fallback_resolution = 1200
)

# Displaying the Farey series values versus the respective index ---------------
# black & white, PDF+EPS
N <- 10L
farey <- c(sort(unique(unlist(lapply(2:N, function(N) (1:(N - 1)) / N)))), 1)
# slope <- diff(range(farey))/(length(farey)-1)
df <- tibble(farey, index = seq_along(farey))

ggplot(df, aes(index, farey)) +
  geom_smooth(method = lm, se = FALSE, colour = "grey60", alpha = 1 / 3) +
  geom_step() +
  # geom_abline(slope = slope, intercept = min(farey)-slope, alpha = 1/3) +
  geom_rug(alpha = 1 / 2) +
  ylab("Farey series")
ggsave("figure5-farey-vs-index.pdf", width = fig_width, height = fig_height)
ggsave("figure5-farey-vs-index.eps",
       width = fig_width, height = fig_height,
       device = cairo_ps, fallback_resolution = 1200
)

# Displaying the objective function values and points evaluated by fminbnd -----
# black & white, PDF+EPS
df <- tibble(
  index = c(
    # percentile
    49, 80, 31, 19, 38, 26, 33, 32, 32, 32, NA, NA,
    49, 80, 31, 35, 39, 37, 37, 36, NA, NA, NA, NA,
    49, 80, 31, 36, 39, 40, 38, 39, NA, NA, NA, NA,
    49, 80, 31, 17, 38, 37, 35, 36, 36, NA, NA, NA,
    49, 80, 31, 54, 46, 47, 40, 44, 43, 42, 43, 43,
    49, 80, 31, 25, 39, 37, 36, 36, 34, 36, NA, NA,
    49, 80, 31, 37, 39, 37, 34, 36, 36, NA, NA, NA,
    # box
    49, 80, 31, 19, 38, 28, 34, 33, 33, NA, NA, NA,
    49, 80, 31, 21, 38, 38, 42, 40, 39, NA, NA, NA,
    49, 80, 31, 33, 39, 43, 37, 36, 38, 37, NA, NA,
    49, 80, 31, 9, 22, 38, 36, 35, 33, 35, 34, NA,
    49, 80, 31, 45, 41, 37, 40, 43, 42, 41, NA, NA,
    49, 80, 31, 19, 38, 37, 35, 33, 36, 35, NA, NA,
    49, 80, 31, 25, 39, 38, 36, 37, 37, 37, NA, NA
  ),
  value = c(
    # percentile
    0.966804, 1.70041, 0.279162, 1.95558, 0.469696, 0.696324,
    0.281147, 0.269359, 0.253947, 0.256774, NA, NA,
    1.16107, 1.61375, 1.11607, 1.03598, 1.03496, 1.0247,
    1.01441, 1.02406, NA, NA, NA, NA,
    1.13347, 1.60767, 1.09458, 0.99964, 0.983566, 1.00411,
    0.989483, 0.984994, NA, NA, NA, NA,
    1.23641, 1.77139, 1.07557, 2.75285, 1.02982, 1.01799,
    1.01941, 1.00784, 1.01448, NA, NA, NA,
    0.980741, 1.13966, 1.1042, 1.00604, 0.976845, 0.978823,
    0.988113, 0.973302, 0.970966, 0.972397, 0.968197, 0.97607,
    1.22222, 1.77281, 1.09174, 1.53933, 1.03222, 1.01024,
    1.0065, 1.01721, 1.01609, 1.01908, NA, NA,
    1.12853, 1.67973, 1.08747, 0.961632, 0.974212, 0.966136,
    1.00251, 0.966865, 0.975122, NA, NA, NA,
    # box
    0.90882, 1.65697, 0.26782, 2.0434, 0.394529, 0.59442,
    0.232481, 0.21761, 0.219033, NA, NA, NA,
    1.02022, 1.5397, 0.881521, 1.90207, 0.790213, 0.776752,
    0.861437, 0.786795, 0.795806, NA, NA, NA,
    0.963359, 1.52588, 0.883782, 0.806283, 0.771955, 0.874007,
    0.724388, 0.759079, 0.764031, 0.718773, NA, NA,
    1.05635, 1.70226, 0.833226, 5.69515, 1.72328, 0.802583,
    0.776241, 0.763807, 0.783029, 0.725408, 0.745387, NA,
    0.835886, 1.06515, 0.867843, 0.798652, 0.728754, 0.754875,
    0.748273, 0.777165, 0.763772, 0.730864, NA, NA,
    1.09691, 1.70285, 0.865888, 2.35372, 0.775974, 0.750855,
    0.725661, 0.807632, 0.787279, 0.78572, NA, NA,
    1.01372, 1.61387, 0.871283, 1.39929, 0.77535, 0.73425,
    0.752021, 0.71723, 0.741082, 0.717375, NA, NA
  ),
  Leaveout = rep(paste("leave-out point", rep(1:7, rep(12, 7))), 2),
  method = rep(c("percentile", "box"), each = 84)
)

ggplot(df, aes(index, value)) +
  geom_point(size = 4, alpha = 1 / 3) +
  xlim(25, 50) +
  ylim(0.2, 1.2) +
  facet_grid(vars(Leaveout), vars(method)) +
  xlab("index of terms in the Farey series of order 20") +
  ylab(expression("value of" ~~ g[N]))
ggsave("figure8-funval-check.pdf", width = fig_width, height = 2.5 * fig_height)
ggsave("figure8-funval-check.eps",
       width = fig_width, height = 2.5 * fig_height,
       device = cairo_ps, fallback_resolution = 1200
)

# Plot model, model+discrepancy function, and experimental values vs. grammage
# colour+black & white, PDF+EPS
model_prctile <- diag(as.matrix(thickness_prctile[, 2 + (1:nobs)]))
model_box <- diag(as.matrix(thickness_box[, 2 + (1:nobs)]))
df <- rbind(
  cbind(
    grammage = thickness_prctile[, "grammage"],
    model = model_prctile,
    "model + discrepancy" = model_prctile + discrepancy_prctile,
    method = "percentile"
  ),
  cbind(
    grammage = thickness_box[, "grammage"],
    model = model_box,
    "model + discrepancy" = model_box + discrepancy_box,
    method = "box"
  ),
  make.row.names = FALSE
)
df <- pivot_longer(
  df,
  cols = starts_with("model"), names_to = "id", values_to = "thickness"
)
df1 <- dplyr::select(thickness, !starts_with("ypred"))
df1 <- pivot_longer(
  df1,
  cols = c(experimental, intrinsic), names_to = "id", values_to = "thickness"
)

p <- ggplot(df, aes(grammage, thickness, colour = id, shape = id)) +
  geom_pointrange(
    aes(grammage, thickness, ymin = thickness - se, ymax = thickness + se),
    data = df1
  ) +
  facet_wrap(vars(method)) +
  geom_point(size = 2) +
  scale_colour_brewer(palette = "Set1") +
  xlab(expression("grammage" ~~ (g / m^2))) +
  ylab(expression("thickness" ~~ (mu * m))) +
  labs(colour = NULL, shape = NULL) +
  guides(
    colour = guide_legend(
      override.aes = list(linetype = c(1, 1, 0, 0))
    )
  ) +
  theme(
    legend.position = "none",
    # legend.position = "top",
    # legend.justification = "right"
  )
p
ggsave("figure9-LOOCV_thickness_unpressed-unrefined_colour.pdf",
  width = fig_width, height = fig_height
)
ggsave("figure9-LOOCV_thickness_unpressed-unrefined_colour.eps",
  width = fig_width, height = fig_height,
  device = cairo_ps, fallback_resolution = 1200
)
p + scale_colour_grey(start = 0, end = 0.5)
ggsave("figure9-LOOCV_thickness_unpressed-unrefined_bw.pdf",
  width = fig_width, height = fig_height
)
ggsave("figure9-LOOCV_thickness_unpressed-unrefined_bw.eps",
  width = fig_width, height = fig_height,
  device = cairo_ps, fallback_resolution = 1200
)

# Displaying the individual discrepancy functions ------------------------------
# colour+black & white, PDF+EPS
n <- 51L
x <- seq(min(thickness$grammage), max(thickness$grammage), length.out = n)
discrepancy_prctile <- as.data.frame(
  vapply(1:nobs, function(ind) SSS_thickness_prctile[[ind]]$Spline(x), double(n))
)
discrepancy_box <- as.data.frame(
  vapply(1:nobs, function(ind) SSS_thickness_box[[ind]]$Spline(x), double(n))
)
df <- as_tibble(
  rbind(
    cbind(x, discrepancy_prctile, method = "percentile"),
    cbind(x, discrepancy_box, method = "box"),
    make.row.names = FALSE
  )
)
names(df)[1:(nobs + 1)] <- c("grammage", paste("point", 1:nobs))
df <- pivot_longer(
  df,
  cols = starts_with("point"), names_to = "Leaveout", values_to = "thickness"
)
df$lab <- ""
df$lab[c(n * nobs - nobs:1, 2 * n * nobs - nobs:1) + 1] <- rep(paste("point", 1:nobs), 2)
df$lab[c(355, 712)] <- "Leave-out point 5"

p <- ggplot(df, aes(grammage, thickness, colour = Leaveout, linetype = Leaveout, label = lab)) +
  geom_line() +
  geom_text_repel() +
  facet_wrap(vars(method)) +
  scale_colour_brewer(palette = "Dark2") +
  geom_hline(yintercept = 0, alpha = 1 / 2) +
  xlab(expression("grammage" ~ ~ (g / m^2))) +
  ylab(expression("thickness" ~ ~ (mu * m))) +
  labs(colour = "Leave-out", linetype = "Leave-out") +
  theme(legend.position = "none")
p
ggsave("figure10-discrepancyfun_thickness_unpressed-unrefined_colour.pdf",
  width = 1.25 * fig_width, height = 1.25 * fig_height
)
ggsave("figure10-discrepancyfun_thickness_unpressed-unrefined_colour.eps",
  width = 1.25 * fig_width, height = 1.25 * fig_height,
  device = cairo_ps, fallback_resolution = 1200
)
p + scale_colour_grey(start = 0, end = 0.5)
ggsave("figure10-discrepancyfun_thickness_unpressed-unrefined_bw.pdf",
  width = 1.25 * fig_width, height = 1.25 * fig_height
)
ggsave("figure10-discrepancyfun_thickness_unpressed-unrefined_bw.eps",
  width = 1.25 * fig_width, height = 1.25 * fig_height,
  device = cairo_ps, fallback_resolution = 1200
)
