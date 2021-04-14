##############################
##############################
#### analyse_example_thermal_niceh.R
#### Edward Lavender (el72@st-andrews.ac.uk)

#### This code:
# 1) Plots and example thermal niche (temperature-abundance curve)
# ... This shows how if the temperature increases,
# ... abundance will go up on cold half of distribution
# ... and down in the warm half of the distribution. 

#### Steps preceding this code:
# 1) Definition of thermal affinities in spptraits dataset (process_spptraits.R)


##############################
##############################
#### Set up 

#### Wipe workspace
rm(list = ls())

#### Load data
spptraits <- readRDS("./data/spptraits.rds")


##############################
##############################
#### Define thermal affinities

#### Pick example species
sp <- "Acipenser oxyrinchus"

#### Get thermal affinities
# Get thermal param 
t10 <- spptraits[spptraits$spp == sp,]$sst_t10
t90 <- spptraits[spptraits$spp == sp,]$sst_t90
sti <- spptraits[spptraits$spp == sp,]$sst_t50
# Define STR
str <- t90 - t10
# Define the STR SD 
str_sd <- (t90-t10)/(2*1.281560031) 


##############################
##############################
#### Graph 

#### Define graphical param
xlim1 <- -10
xlim2 <- 40

#### Set up to save figure
tiff("./fig/example_thermal_niche.tiff", 
     height = 4, width = 5, units = "in", res = 600)

#### Make blank plot of abundance ~ temperature 
# 'Abundance' will be scaled later so we'll set the axes to be 0 - 1 here 
plot(NA,
     axes = FALSE,
     xlab = "", ylab = "",
     xlim = c(xlim1, xlim2), ylim = c(0, 1), 
     type = "n")

#### Add scaled abundance
# Define sequence of temperatures and calculate scaled abundance based on thermal parameters 
x <- seq(xlim1, xlim2, by = 0.1)
curve(dnorm(x, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd), 
      add = TRUE, lty = 1, lwd = 0.5)

#### Add sti
lines(c(sti, sti), c(0, 1), lty = 2)
px <- par(xpd = NA)
text(sti, 1+0.05, "STI")
par(px)

#### Colour the warm/cold halfs of the distribution 
position1 <- which(x < sti)[length(which(x < sti))]
position2 <- which(x > sti)[1]
x_cold <- x[1:position1]
x_warm <- x[position2:length(x)]
y <- dnorm(x, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd)
y_cold <- y[1:position1]
y_warm <- y[position2:length(y)]
polygon(c(x_cold, rev(x_cold)), c(y_cold, rep(0, length(y_cold))), 
        col = scales::alpha("skyblue", 0.3), border = NA)
polygon(c(x_warm, rev(x_warm)), c(y_warm, rep(0, length(y_warm))), 
        col = scales::alpha("red", 0.3), border = NA)

#### Example changes in warm half of distribution
# Define two temperatures in warm half, before and after climate change 
t1 <- round(sti + 5, digits = 0); t1 # before
t2 <- round(t1 + 2, digits = 0); t2  # after hypothetical climate change 
# Calculate the abundance at each temperature
abundance_t1 <- dnorm(t1, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd)
abundance_t2 <- dnorm(t2, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd)
# Add the abundance at temp1 and temp2 to the graph
points(t1, abundance_t1, pch = 23, bg = "black")
points(t2, abundance_t2, pch = 21, bg = "black")
# Add arrows defining the direction of change 
offset <- 0.05
lines(c(t1, t2 + offset), c(abundance_t1, abundance_t1))
arrows(x0 = t2 + offset, x1 = t2 + offset, y0 = abundance_t1, y1 = abundance_t2 + 0.025, length = 0.05)

#### Example changes in cold half of distribution 
# Define two temperatures in cold half, before and after climate change 
t1 <- round(sti - 5, digits = 0); t1  # before
t2 <- round(t1 + 2, digits = 0); t2   # after climate change
# Calculate abundance at each temperature 
abundance_t1 <- dnorm(t1, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd)
abundance_t2 <- dnorm(t2, mean = sti, sd = str_sd)/dnorm(sti, mean = sti, sd = str_sd)
# Add to plot 
points(t1, abundance_t1, pch = 23, bg = "black")
points(t2, abundance_t2, pch = 21, bg = "black")
# Add arrows showing the direction of change 
lines(c(t1, t2 + offset), c(abundance_t1, abundance_t1))
arrows(x0 = t2 + offset, x1 = t2 + offset, y0 = abundance_t1, y1 = abundance_t2 - 0.025, length = 0.05)

#### Add legend
legend("topright", pch = c(23, 21), pt.bg = c("black", "black"),  col = c("black", "black"), 
       legend = c(expression(IRA[T[0]]), 
                  expression(IRA[T[0] ~ paste("+", Delta, "T", sep = "")])), 
       box.lty = 3,
       box.lwd = 0.5,
       y.intersp = 1.5)

#### Add labels 
axis(side = 1, seq(xlim1, xlim2, by = 10), pos = 0)
axis(side = 2, seq(0, 1, by = 0.25), pos = xlim1, las = 2)
mtext(side = 1, expression(paste("Temperature, T (", degree, "C)", sep = "")), line = 2.5)
mtext(side = 2, "Relative Abundance, A", line = 2.5)

#### Save 
dev.off()


#### End of code. 
##############################
##############################