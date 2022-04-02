library(tibble)
library(dplyr)
library(ggplot2)
library(readr)

# Experimental data ------------------------------------------------------------
# MorFi Lab unit
coarseness <- mean(c(0.0938, 0.0965, 0.0964)) # mg/m
fibre_width <- mean(c(18.8, 18.8, 18.8)) # um
fibre_length <- mean(c(797, 801, 799)) # um
# from SEM images, um
fibre_thickness <- c(
  7.356, 9.255, 5.774, 5.715, 6.882, 7.360, 7.356, 5.700, 3.797, 5.932,
  8.546, 10.939, 9.267, 8.546, 3.559, 5.735, 8.546, 4.034, 9.255, 4.746,
  5.989, 5.519, 4.041, 2.556, 5.253, 8.072, 5.715, 5.422, 3.559, 4.458,
  4.989, 4.515, 4.534, 3.752, 4.662, 5.640, 6.050, 4.363, 2.848, 6.124,
  6.211, 5.306, 5.210, 4.477, 5.695, 4.041, 2.973, 4.805, 5.463, 4.376,
  2.373, 1.186, 2.857, 3.504, 3.373, 1.898, 4.041, 4.515, 2.443, 3.752,
  4.680, 4.245, 5.221, 4.245, 4.502, 3.166, 3.422, 3.834, 3.297, 2.857,
  3.094, 3.630, 5.870, 4.278, 3.020, 3.834, 4.096, 4.144, 3.949, 4.509,
  5.269, 3.039, 4.198, 4.805, 4.278, 5.411, 5.242, 3.707, 1.384, 3.192,
  6.649, 2.522, 2.706, 4.892, 3.614, 2.706, 3.305, 3.892, 3.528, 2.027,
  1.678, 2.567, 1.853, 3.630, 3.297, 4.083, 5.540, 3.559, 4.414, 4.041,
  4.989, 3.622, 2.251, 2.041, 2.706, 5.793, 4.775, 2.857, 4.869, 3.039,
  3.949, 4.041, 2.767, 3.528, 2.935, 4.458, 3.373, 6.188, 2.443, 2.337,
  6.064, 2.477, 1.278, 5.715, 3.684, 3.978, 3.002, 2.556, 3.085, 3.691,
  2.706, 6.050, 5.006, 2.443, 4.414, 2.653, 2.621, 2.149, 1.853, 2.477,
  4.041, 3.863, 3.228, 1.913, 3.020, 3.892, 4.583, 2.567, 3.002, 3.892,
  1.424, 6.979, 2.522, 2.777, 6.050, 2.897, 2.621, 5.204, 3.559, 2.706,
  3.622, 3.528, 4.278, 2.014, 7.582, 4.775, 2.935, 4.414, 3.707, 2.522,
  3.949, 3.614, 2.239, 1.210, 3.085, 3.622, 3.356, 4.198, 3.707, 2.706,
  3.528, 1.898, 3.528, 4.955, 4.151, 4.515, 2.653, 4.989, 5.034, 4.607,
  3.978, 3.002, 4.041, 4.198, 3.804, 3.567, 5.416, 2.857, 3.121, 5.269,
  6.563, 1.853, 2.621, 4.271, 3.121, 1.210, 2.973, 1.913, 1.913, 4.607,
  3.559, 2.897, 2.935, 2.394, 4.938, 1.061, 4.278, 3.630, 3.559, 2.897,
  4.083, 6.425, 2.767, 4.034, 5.937, 1.501, 1.913, 3.949, 2.848, 4.746,
  8.595, 2.027, 5.380, 2.443, 2.443, 3.707, 3.398, 2.041, 3.192, 4.515,
  4.886, 3.305, 2.027, 2.027, 3.978, 1.913, 3.002, 2.621, 2.373, 2.251,
  1.186, 4.477, 2.767, 5.463, 4.583, 5.837, 5.354, 3.504, 2.556, 2.014,
  2.477, 5.028, 4.892, 4.509, 2.767, 4.245, 3.039, 6.211, 3.094, 2.867,
  3.914, 2.349, 4.564, 3.121, 3.630, 4.886, 6.318, 4.205, 5.715, 3.094,
  3.826, 4.034, 3.305, 3.455, 1.443, 4.151, 2.973, 3.166, 2.522, 3.622,
  2.556, 3.455, 2.653, 3.504, 4.477, 4.746, 2.857, 2.522, 6.367, 3.892,
  3.559, 1.661, 1.711, 3.192, 3.085, 5.932, 2.767, 2.027, 3.826, 1.913,
  2.777, 5.779
)
density_fibre <- coarseness / (fibre_width * mean(fibre_thickness)) * 1e3 # g/cm3
sheet_former_area <- 0.02138 # m2
nobs <- 7
repeats <- 10

# unrefined pulp, no pressing
handsheet <- tibble(
  micrometer = rep("FRANK-PTI", each = repeats * nobs),
  pressing = rep("no pressing", each = repeats * nobs),
  refining = rep("no", each = repeats * nobs),
  id = rep(1:nobs, each = repeats),
  # Mettler Toledo - PB303 DeltaRange, g
  mass = c(
    0.437, 0.444, 0.456, 0.458, 0.468, 0.452, 0.441, 0.461, 0.465, 0.437,
    0.911, 0.911, 0.911, 0.912, 0.902, 0.907, 0.883, 1.022, 0.740, 0.968,
    1.354, 1.369, 1.337, 1.364, 1.370, 1.378, 1.368, 1.376, 1.381, 1.374,
    1.862, 1.864, 1.870, 1.878, 1.925, 1.880, 1.858, 1.802, 1.841, 1.830,
    2.323, 2.343, 2.350, 2.344, 2.289, 2.350, 2.339, 2.308, 2.339, 2.347,
    2.747, 2.764, 2.745, 2.725, 2.730, 2.774, 2.748, 2.756, 2.817, 2.773,
    3.254, 3.238, 3.242, 3.252, 3.187, 3.172, 3.227, 3.242, 3.227, 3.218
  ),
  # FRANK-PTI, tissue acc. to ISO 12625-3, um
  thickness = c(
    109, 104, 113, 104, 113, 112, 107, 110, 115, 111,
    175, 174, 185, 185, 176, 181, 168, 191, 157, 195,
    269, 265, 257, 258, 247, 238, 233, 239, 228, 235,
    323, 309, 326, 341, 315, 328, 312, 313, 323, 322,
    422, 425, 408, 409, 419, 405, 416, 427, 427, 416,
    492, 459, 474, 512, 496, 475, 484, 474, 497, 494,
    613, 544, 573, 540, 593, 517, 566, 549, 566, 568
  ),
  grammage = mass / sheet_former_area, # g/m2
  porosity = 1 - grammage / (thickness * density_fibre)
)
# Grammage as the aggregate from individual sheets
handsheet <- mutate(
  handsheet,
  grammage_avg = rep(
    pull(summarise(group_by(handsheet, id), mean(grammage))),
    each = repeats
  )
)

# Indirect method used for obtaining the intrinsic thickness and the intrinsic porosity
w <- rep(
  pull(summarise(group_by(handsheet, id), sd(thickness))),
  each = repeats
)
w <- 1/w^2
correct <- lm(thickness ~ grammage_avg, data = handsheet, weights = w)

# Plot settings ----------------------------------------------------------------
fig_width <- 7
fig_asp <- .75
fig_height <- fig_asp * fig_width

# Plots for the the calibration procedure --------------------------------------
# Empirical density of fibre thickness - black & white, PDF+EPS
ggplot(enframe(fibre_thickness, name = NULL, value = "fibre_thickness"), aes(fibre_thickness)) +
  geom_density(bw = "SJ") +
  geom_vline(xintercept = mean(fibre_thickness), alpha = 1 / 2) +
  xlab(expression("fibre thickness" ~~ (mu * m)))
ggsave("figure6-fibre-thickness-expt_unpressed-unrefined.pdf",
  width = fig_width, height = fig_height
)
ggsave("figure6-fibre-thickness-expt_unpressed-unrefined.eps",
  width = fig_width, height = fig_height,
  device = cairo_ps, fallback_resolution = 1200
)

# Sheet thickness vs. grammage - colour+black & white, PDF+EPS
df <- tibble(
  grammage_avg = pull(summarise(group_by(handsheet, id), mean(grammage))),
  thickness_avg = pull(summarise(group_by(handsheet, id), mean(thickness)))
)
ggplot(handsheet, aes(grammage_avg, thickness)) +
  geom_point() +
  geom_smooth(method = lm) +
  # stat_summary(geom = "point", fun = "mean", size = 4, alpha = 1 / 3) +
  geom_point(aes(grammage_avg, thickness_avg), data = df, size = 4, alpha = 1 / 3) +
  xlab(expression("grammage" ~~ (g / m^2))) +
  ylab(expression("thickness" ~~ (mu * m)))
ggsave("figure7-sheet-thickness-expt-tissue_unpressed-unrefined_colour.pdf",
  width = fig_width, height = fig_height
)
ggsave("figure7-sheet-thickness-expt-tissue_unpressed-unrefined_colour.eps",
  width = fig_width, height = fig_height,
  device = cairo_ps, fallback_resolution = 1200
)
ggplot(handsheet, aes(grammage_avg, thickness)) +
  geom_point() +
  geom_smooth(colour = "grey40", method = lm) +
  # stat_summary(geom = "point", fun = "mean", size = 4, alpha = 1 / 3) +
  geom_point(aes(grammage_avg, thickness_avg), data = df, size = 4, alpha = 1 / 3) +
  xlab(expression("grammage" ~~ (g / m^2))) +
  ylab(expression("thickness" ~~ (mu * m)))
ggsave("figure7-sheet-thickness-expt-tissue_unpressed-unrefined_bw.pdf",
  width = fig_width, height = fig_height
)
ggsave("figure7-sheet-thickness-expt-tissue_unpressed-unrefined_bw.eps",
  width = fig_width, height = fig_height,
  device = cairo_ps, fallback_resolution = 1200
)

# Data transfer for the the calibration procedure ------------------------------
write_csv(tibble(
  thickness = mean(fibre_thickness),
  width = fibre_width,
  length = fibre_length,
  coarseness = coarseness
), "fibre_morphology.dat")

# -> MATLAB script 'calibration_suppinfo.m'
# Experimental data
exptdata2cal <- summarise(
  group_by(handsheet, id),
  grammage_avg = mean(grammage), grammage_stddev = sd(grammage),
  thickness_avg = mean(thickness), thickness_stddev = sd(thickness),
  porosity_avg = mean(porosity), porosity_stddev = sd(porosity)
)
exptdata2cal <- mutate(
  exptdata2cal,
  thickness_correct_avg = thickness_avg - coef(correct)[["(Intercept)"]],
  porosity_correct = 1 - 1/(coef(correct)[[2]] * density_fibre)
)

write_csv(exptdata2cal, "exptdata2cal_thickness-tissue-ISO-12625-3_grammage-all_unpressed-unrefined.dat")
