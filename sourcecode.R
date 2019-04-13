## ------------------------------------------------------------------------------------------------#
## Before running the script, please set the working directory to the script directory
## ------------------------------------------------------------------------------------------------#

library("dMod")
library("cowplot")
library("parallel")

mywd <- getwd()

if (!dir.exists("Computations")) dir.create("Computations")
setwd("Computations")

# 4.1 Simulation and prediction  -------------------------------------------------------------------


reactions <- NULL
reactions <- addReaction(reactions, "TCA_buffer", "TCA_cell",
                         rate = "import*TCA_buffer",
                         description = "Uptake")
reactions <- addReaction(reactions, "TCA_cell", "TCA_buffer",
                         rate = "export_sinus*TCA_cell",
                         description = "Sinusoidal export")
reactions <- addReaction(reactions, "TCA_cell", "TCA_cana",
                         rate = "export_cana*TCA_cell",
                         description = "Canalicular export")
reactions <- addReaction(reactions, "TCA_cana", "TCA_buffer",
                         rate = "reflux*TCA_cana",
                         description = "Reflux into the buffer")

mymodel <- odemodel(reactions, modelname = "bamodel")

x <- Xs(mymodel, condition = NULL)

times <- seq(0, 50, .1)
pars <- c(TCA_buffer = 1,
          TCA_cell = 0,
          TCA_cana = 0,
          import = 0.2,
          export_sinus = 0.2,
          export_cana = 0.04,
          reflux = 0.1)

out <- x(times, pars)
plot(out)

outSens <- getDerivs(x(times, pars[4:7], fixed = pars[1:3]))
plot(outSens)

# Reproduce Figure 2 in the paper:
P.pred <- plot(out) + 
  facet_wrap(~name, scales="free", ncol=1) + 
  theme(legend.position = "none") + ylab("concentration")
P.sens <- plot(outSens) + 
  facet_wrap(~name, scales="free", ncol=3) + 
  theme(legend.position = "none") + ylab("parameter sensitivity")

P <- plot_grid(P.pred, P.sens, 
               rel_widths = c(0.25, 0.75),
               labels = c("A", "B"), label_size = 10)

ggsave("../Figures/figure1.pdf", P, 
       height=10.9/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")


# 4.2 Observation function and simulated data ------------------------------------------------------


observables <- eqnvec(
  buffer = "s*TCA_buffer",
  cellular = "s*(TCA_cana + TCA_cell)"
)
g <- Y(observables, f = x, condition = NULL,
       compile = TRUE, modelname = "obsfn")

pars["TCA_cell"] <- 0.3846154
pars["TCA_cana"] <- 0.1538462
pars["TCA_buffer"] <- 0
pars["s"] <- 1e3

out <- (g*x)(times, pars, conditions = "standard")

suppressWarnings(RNGversion("3.5.0"))
set.seed(1)

timesD <- c(0.1, 1, 3, 7, 11, 15, 20, 41)
datasheet <- subset(as.data.frame(out),
                    time %in% timesD & name %in% names(observables))

datasheet <- within(datasheet, {
  sigma <- sqrt(value + 1)
  value <- rnorm(length(value), value, sigma)
})

data <- as.datalist(datasheet)
plot(out, data)

# Reproduce Figure 3 in the paper:
P.simu <- plot(out, data) + 
  facet_wrap(~name,scales="free", nrow=1) + 
  theme(legend.position = "none") + 
  ylab("concentration")

ggsave("../Figures/figure2.pdf", P.simu, 
       height=5.1/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")


# 4.3 Parameter transformation ---------------------------------------------------------------------

p <- P(
  trafo = eqnvec(
    TCA_buffer = "0",
    TCA_cell = "exp(TCA_cell)",
    TCA_cana = "exp(TCA_cana)",
    import = "exp(import)",
    export_sinus = "exp(export_sinus)",
    export_cana = "exp(export_cana)",
    reflux = "exp(reflux)",
    s = "exp(s)"
  ),
  condition = "standard"
)
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
plot((g*x*p)(times, pouter), data)

# Reproduce Figure 4 in the paper
P.fit <- plot((g*x*p)(times, pouter), data) + 
  facet_wrap(~name,scales="free", nrow=1) + 
  theme(legend.position = "none") + 
  ylab("concentration")


ggsave("../Figures/figure3.pdf", P.fit, 
       height=5.1/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# 4.4 Objective function and model fitting ---------------------------------------------------------

obj <- normL2(data, g*x*p) + constraintL2(pouter, sigma = 10)
myfit <- trust(obj, pouter, rinit = 1, rmax = 10)

plot((g*x*p)(times, myfit$argument), data)

# Reproduce Figure 5 in the paper
P.fit <- plot((g*x*p)(times, myfit$argument), data) + 
  facet_wrap(~name,scales="free", nrow=1) + 
  theme(legend.position = "none") + 
  ylab("concentration")
ggsave("../Figures/figure4.pdf", P.fit,
       height=5.1/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")


out_mstrust <- mstrust(obj, pouter, sd = 4, fits = 50,
                       rinit = 1, rmax = 10, iterlim = 500,
                       cores = 4)

myframe <- as.parframe(out_mstrust)

plotValues(myframe, tol = .01, value < 100)
plotPars(myframe, tol = .01, value < 100)

# Reproduce Figure 6 in the paper
P.lhs <- plotValues(myframe, tol = 0.01, value < 100) + 
  ylab(expression(chi^2~value)) 

P.pars <- plotPars(myframe, tol = 0.01,value < 100) + 
  xlab("")  + 
  ylab("parameter value")

uniquepars <- unique(subset(myframe, value < 100), tol = 0.01)
myprediction <- predict(g*x*p, pars = uniquepars, times = times)
P.pred <- ggplot() + 
  geom_line(data=myprediction, aes(x=time, y=value, color=factor(.value))) + 
  geom_pointrange(data=data[[1]], 
                  aes(x=time, y=value, ymin=value - sigma, ymax=value + sigma), shape=20) +
  facet_wrap(~name, scales="free", nrow=1) + theme_dMod() + theme(legend.position = "none") + 
  scale_color_dMod() + ylab("concentration")


P <- ggdraw() + 
  draw_plot(P.lhs, 0, .4, .5, .6) +
  draw_plot(P.pars, .5, .4, .5, .6) +
  draw_plot(P.pred, 0, 0, 1, .4) +
  draw_plot_label(c("A", "B", "C"), c(0,.5, 0), c(1,1,.4), size = 10)

ggsave("../Figures/figure5.pdf", P, height=12.7/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# 4.5 Working with several conditions --------------------------------------------------------------

pars["reflux"] <- 1e3
out <- (g*x)(times, pars, conditions = "open")
datasheet <- subset(as.data.frame(out),
                    time %in% timesD & name %in% names(observables))
datasheet <- within(datasheet, {
  sigma <- sqrt(value + 1)
  value <- rnorm(length(value), value, sigma)
})

data <- data + as.datalist(datasheet)

trafo <- getEquations(p, conditions = "standard")
trafo["reflux"] <- "exp(reflux_open)"
p <- p + P(trafo, condition = "open")

outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

obj <- normL2(data, g*x*p) + constraintL2(pouter, sigma = 10)

out_mstrust <- mstrust(obj, pouter, rinit = 1, rmax = 10, iterlim = 500,
                       sd = 4, cores = 4, fits = 50)
myframe <- as.parframe(out_mstrust)
plotValues(myframe, tol = 0.1, value < 100)
plotPars(myframe, tol = 0.1, value < 100)

bestfit <- as.parvec(myframe)
plot((g*x*p)(times, bestfit), data)

# Reproduce Figure 7 in the paper
d.pred <- as.data.frame((g*x*p)(times, bestfit))
d.data <- as.data.frame(data)
P.fit <- 
  ggplot(mapping = aes(x = time, y = value, 
                       ymin = value - sigma, ymax = value + sigma)) +
  facet_wrap(~name, scales = "free", nrow = 1) +
  geom_errorbar(data = d.data, width = 0) +
  geom_line(data = d.pred, aes(lty = condition, pch = condition)) +
  geom_point(data = d.data, aes(lty = condition, pch = condition)) +
  scale_shape_manual(name = "condition",
                     values = c(standard = 19, open = 1)) +
  scale_linetype_manual(name = "condition",
                        values = c(standard = 1, open = 2)) +
  ylab("concentration") +
  theme_dMod() +
  theme(legend.key.width = unit(1, "cm"))
  
P.lhs <- plotValues(myframe, tol = 0.01, value < 100) +
  ylab(expression(chi^2~value))

P.pars <- plotPars(myframe, tol = 0.01,value < 100) + 
  xlab("")  + ylab("parameter value")


P <- ggdraw() + 
  draw_plot(P.fit, 0, .6, 1, .4) +
  draw_plot(P.lhs, 0, 0, .5, .6) +
  draw_plot(P.pars, .5, 0, .5, .6) +
  draw_plot_label(c("A", "B", "C"), c(0,0, .5), c(1,.6,.6), size = 10)

ggsave("../Figures/figure6.pdf", P, height=11.4/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# 4.6 Parameter uncertainty and identifiability ----------------------------------------------------

profiles <- profile(obj, bestfit, names(bestfit),
                    limits = c(-5, 5), cores = 4)
plotProfile(profiles)
plotPaths(profiles, whichPar = "s")

# Reproduce Figure 8 in the paper
P.prof <- plotProfile(profiles) + 
  facet_wrap(~name, scales="free_x", nrow=2) + 
  guides(color = "none") +
  theme(legend.position = "top")

P.path <- plotPaths(profiles, whichPar = "s") + 
  facet_wrap(~combination, nrow=1) + 
  theme(legend.position = "none")

P <- ggdraw() + 
  draw_plot(P.prof, 0, .3, 1, .7) +
  draw_plot(P.path, 0, 0, 1, .3) +
  draw_plot_label(c("A", "B"), c(0,0), c(1,.3), size = 10)

ggsave("../Figures/figure7.pdf", P, height=15.2/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# Save profiles of model without steady state assumption for the next section
profiles_noSS <- profiles


# 4.7 Steady-state constraints and implicit transformations ----------------------------------------

pSS <- NULL
equations <- getEquations(p)
conditions <- names(equations)
for (n in conditions) {
  equations[[n]]["TCA_cana"] <- "exp(export_cana)*exp(TCA_cell)/exp(reflux)"
  pSS <- pSS + P(equations[[n]], condition = n)
}

# Produce profiles with explicit steady state parameterization (code not shown in paper)
outerpars <- getParameters(pSS)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)

obj <- normL2(data, g*x*pSS) + constraintL2(pouter, sigma = 10)
bestfit <- trust(obj, pouter, rinit = 1, rmax = 10)$argument

profiles_SS_analytic <- profile(obj, bestfit, names(bestfit), limits = c(-5, 5), cores = 4)


# Implementation of the system with implicit steady-state transformation
# Construct reactions
reactions <- NULL
reactions <- addReaction(reactions, "TCA_buffer", "TCA_cell",
                         rate = "import*TCA_buffer",
                         description = "Uptake")
reactions <- addReaction(reactions, "TCA_cell", "TCA_buffer",
                         rate = "export_sinus*TCA_cell",
                         description = "Sinusoidal export")
reactions <- addReaction(reactions, "TCA_cell", "TCA_cana",
                         rate = "export_cana*TCA_cell",
                         description = "Canalicular export")
reactions <- addReaction(reactions, "TCA_cana", "TCA_buffer",
                         rate = "(reflux*(1-switch) + reflux_open*switch)*TCA_cana",
                         description = "Reflux into the buffer")
reactions <- addReaction(reactions, "0", "switch",
                         rate = "0",
                         description = "Create a switch")

# Construct events
events <- NULL
events <- addEvent(events, var = "TCA_buffer", time = 0, value = 0      )
events <- addEvent(events, var = "switch"    , time = 0, value = "OnOff")

# Translate into ODE model
mymodel <- odemodel(reactions, modelname = "bamodel2", events = events)

# Set up prediction function
x <- Xs(mymodel)


# Set up implicit parameter transformation
f <- as.eqnvec(reactions)[c("TCA_buffer", "TCA_cana", "TCA_cell")]
f["TCA_cell"] <- "TCA_buffer + TCA_cana + TCA_cell - TCA_tot"
pSS <- P(f, method = "implicit",
         compile = TRUE, modelname = "pfn")

# Set up explicit parameter transformation
innerpars <- unique(c(getParameters(mymodel),
                      getSymbols(observables),
                      getSymbols(f)))

# Use repar() to replace left-hand side by right-hand side of formula
trafo <- repar("x~x"     , x = innerpars)
trafo <- repar("x~0"     , x = reactions$states, trafo)
trafo <- repar("x~exp(x)", x = setdiff(innerpars, "OnOff"), trafo)

p <- 
  P(repar("OnOff~0", trafo), condition = "standard") +
  P(repar("OnOff~1", trafo), condition = "open")
  

# Generate observation function with modified equations
g <- Y(observables, f = x,
       compile = TRUE, modelname = "obsfn2")

# Generate objective function
outerpars <- getParameters(p)
pouter <- structure(rep(-1, length(outerpars)), names = outerpars)
obj <- normL2(data, g*x*pSS*p) + constraintL2(pouter, sigma = 10)

# Search for best fit and compute profiles (code not shown in paper)
out_mstrust <- mstrust(obj, pouter, rinit = 1, rmax = 10, iterlim = 500,
                       sd = 4, cores = 4, fits = 50)
myframe <- as.parframe(out_mstrust)
bestfit <- as.parvec(myframe)

profiles <- profile(obj, bestfit, names(bestfit),
                    limits = c(-5, 5), cores = 4)

profiles_SS_implicit <- profiles



# Reproduce Figure 9 in the paper
profiles_tot <- list("noSS" = profiles_noSS, 
                     "SS_explicit" = profiles_SS_analytic,
                     "SS_implicit" = profiles_SS_implicit)

myplot <- plotProfile(profiles_tot, mode == "data")
plotData <- attr(myplot, "data")
plotDataSparse <- do.call(
  rbind, 
  lapply(split(plotData, list(plotData$name, plotData$proflist), drop = TRUE), function(d) {
    N <- 12
    indx.sparse <- seq(1, nrow(d), length.out = N+2)[-c(1, N+2)]
    d[indx.sparse,]
  })
)



P <- ggplot(plotData, aes(x = par, 
                          y = delta, 
                          group = interaction(proflist, mode), 
                          color = proflist, 
                          pch = proflist)) + 
  facet_wrap(~name, scales="free_x", nrow = 2) + 
  geom_hline(yintercept=c(1, 2.71, 3.84), lty=2, color="gray") + 
  geom_line() +
  geom_point(data = subset(plotData, is.zero), pch = 19, show.legend = FALSE) +
  ylab(expression(paste("CL /", Delta*chi^2))) +
  scale_y_continuous(breaks=c(1, 2.7, 3.84), labels = c("68% / 1   ", "90% / 2.71", "95% / 3.84"), 
                     limits = c(NA, 4)) +
  xlab("parameter value") +
  scale_shape_manual(name = "implementation", 
                     values = c(noSS = NA, SS_explicit = 3, SS_implicit = 4),
                     labels = c(noSS = "no SS", 
                                SS_explicit = "explicit SS", 
                                SS_implicit = "implicit SS")) +
  scale_color_dMod(name="implementation",
                   labels = c(noSS = "no SS", 
                              SS_explicit = "explicit SS",
                              SS_implicit = "implicit SS")) + 
  geom_point(data = plotDataSparse) +
  theme_dMod() +
  theme(legend.position = c(0.9, 0.2))



ggsave("../Figures/figure8.pdf", P, 
       height=10.2/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# 4.8 Prediction uncertainty and validation profiles -----------------------------------------------


obj.validation <- datapointL2(name = "TCA_cell",
                              time = 41,
                              value = "d1",
                              sigma = .002,
                              condition = "standard")


fixed <- c(TCA_tot = log(1))
myfit <- trust(obj + obj.validation,
               parinit = c(d1 = 1, pouter[!names(pouter) %in% names(fixed)]),
               fixed = fixed,
               rinit = 1, rmax = 10)

pred <- (g*x*pSS*p)(times = exp(seq(log(1e-1), log(50), length.out = 100)),
                    pars = myfit$argument, 
                    fixed = c(TCA_tot = log(1)))



profile_prediction <- profile(obj + obj.validation,
                              myfit$argument, "d1", limits = c(-5, 5),
                              fixed = fixed)


# Compute validation band to reproduce the Figure 10 in the paper
tPred <-  c(1, 3, 6, 10, 15, 21, 28, 36, 45)
obj <- normL2(data, g*x*pSS*p, times = tPred) + constraintL2(pouter, sigma = 10)

prediction_band <- do.call(rbind, lapply(tPred, function(t) {
  
  obsname <- "TCA_cell"
  condition <- "standard"
  fixed <- c(TCA_tot = log(1))
  parinit <- unclass(bestfit)[setdiff(names(bestfit), names(fixed))]
  
  cat("\rt = ", t)
  
  obj.validation <- datapointL2(name = obsname,
                                time = t,
                                value = "d1",
                                sigma = .002,
                                condition = condition)
  
  
  myfit <- trust(obj + obj.validation,
                 parinit = c(d1 = 1, parinit),
                 fixed = fixed,
                 rinit = 1, rmax = 10)
  
  
  profile_prediction <- suppressMessages(
    profile(obj + obj.validation,
            myfit$argument, "d1", limits = c(-5, 5),
            fixed = c(TCA_tot = log(1)))
  )
  
  conf <- confint(profile_prediction, "d1")
  
  with(conf, {
    data.frame(
      name = obsname,
      time = t,
      value = value,
      lower = lower,
      upper = upper,
      condition = condition
    )
  })
  
}))



P.prof <- plotProfile(profile_prediction, !mode=="prediction") + theme(legend.position = "bottom") + 
  scale_linetype(guide=guide_legend(ncol=1, title.position = "top")) + 
  scale_color_dMod(guide=FALSE) +
  scale_x_continuous(breaks = c(0.19,0.2,0.21)) +
  theme(legend.background = element_rect(fill = "white", color = "black"))


P.fit <- plot(pred, data) +
  scale_x_log10() +
  geom_ribbon(data = prediction_band, aes(ymin = lower, ymax = upper), 
              lty = 0, alpha = .3, show.legend = FALSE) +
  geom_point(data = prediction_band, aes(x = time, y = lower, color = condition), 
             inherit.aes = FALSE, pch = 4, show.legend = FALSE) +
  geom_point(data = prediction_band, aes(x = time, y = upper, color = condition), 
             inherit.aes = FALSE, pch = 4, show.legend = FALSE) +
  geom_line(data = prediction_band, aes(x = time, y = lower, color = condition), 
            inherit.aes = FALSE, lty = 2, show.legend = FALSE) +
  geom_line(data = prediction_band, aes(x = time, y = upper, color = condition), 
            inherit.aes = FALSE, lty = 2, show.legend = FALSE) +
  theme(legend.position = c(0.88, 0.2), 
        legend.background = element_rect(fill = "white", color = "black"))


P <- ggdraw() + 
   draw_plot(P.prof, 0, .1, .3, .9) +
   draw_plot(P.fit, .3, 0, .7, 1) +
   draw_plot_label(c("A", "B"), c(0,.3), c(1,1), size = 10)
  
ggsave("../Figures/figure9.pdf", P, 
       height=10.2/2.54, width=20.3/2.54, 
       device = cairo_pdf, family = "LM Roman 10")

# 4.9 Speed comparison (Code not shown in paper) ---------------------------------------------------

plist <- lapply(1:50, function(i) pouter + rnorm(length(pouter)))
time1 <- mclapply(plist, function(pars) {
  system.time(trust(obj, pars, rinit = 1, rmax = 10))[1]
}, mc.cores = 4, mc.preschedule = FALSE)
time2 <- mclapply(plist, function(pars) {
  system.time(optim(pars,
                    function(...) obj(..., deriv = FALSE)$value,
                    control = list(maxit = 10000)))[1]
}, mc.cores = 4, mc.preschedule = FALSE)
time3 <- mclapply(plist, function(pars) {
  system.time(optim(pars,
                    function(...) obj(..., deriv = FALSE)$value,
                    method = "L-BFGS-B",
                    lower = -10, upper = 10))[1]
}, mc.cores = 4, mc.preschedule = FALSE)

mytable <- rbind(
  data.frame(
    code = factor("compiled", levels = c("compiled", "not compiled")),
    optimizer = factor(c(rep("trust-region", length(time1)),
                         rep("Nelder-Mead", length(time2)),
                         rep("L-BFGS-B", length(time3))),
                       levels = c("trust-region", "Nelder-Mead", "L-BFGS-B")),
    time = c(unlist(time1), unlist(time2), unlist(time3))),
  data.frame(
    code = factor("not compiled", levels = c("compiled", "not compiled")),
    optimizer = factor(c(rep("trust-region", length(time1)),
                         rep("Nelder-Mead", length(time2)),
                         rep("L-BFGS-B", length(time3))),
                       levels = c("trust-region", "Nelder-Mead", "L-BFGS-B")),
    time = 50*c(unlist(time1), unlist(time2), unlist(time3)))
)


P <- ggplot(mytable, aes(x  = optimizer, y = time, color = code, pch = code)) +
  geom_violin(position = "identity", size = .3) + geom_jitter(size = 1, alpha = .3) +
  scale_y_log10(breaks = c(1, 5, 10, 50, 100, 500, 1000, 5000),
                minor_breaks = c(seq(1e-1, 1e0, 1e-1),
                                 seq(1e0 , 1e1, 1e0 ),
                                 seq(1e1 , 1e2, 1e1 ),
                                 seq(1e2 , 1e3, 1e2 ),
                                 seq(1e3 , 1e4, 1e3 ))) +
  scale_x_discrete(labels = c("(1) trust-region", "(2) Nelder-Mead", "(3) L-BFGS-B")) +
  ylab("runtime for a single fit [s] (logarithmic)") +
  xlab("optimization method") +
  scale_shape_manual(values = c(1, 2)) +
  scale_color_manual(values = c("black", "darkgoldenrod3")) +
  annotate(
    geom = "text",
    x = c("trust-region", "Nelder-Mead", "L-BFGS-B", "trust-region", "Nelder-Mead", "L-BFGS-B"),
    y = c(.2, 1.5, 300, 7, 60, 5),
    label = c("(dMod)", "(mkin, FME)*", "(FME)*", "(nlmeODE)*", "(scaRabee, FME)*", "(FME)*"),
    hjust = 0, size = 2, family = "LM Roman 10"
  ) +
  annotate(
    geom = "text", y = .2, x = "L-BFGS-B", label = "* estimated",
    hjust = 0, vjust = -1, size = 2, family = "LM Roman 10", fontface = "bold"
  ) +
  theme_dMod(base_size = 8, base_family = "LM Roman 10") +
  theme(legend.position = "top",
        legend.background = element_blank(),
        legend.box.margin=margin(-5,-5,-5,-5),
        axis.text.y = element_text(hjust = 0)) +
  coord_flip()


cairo_pdf(filename = "../Figures/figure10.pdf",
          width = 15.2/2.54, height = 5.1/2.54)
print(P)
dev.off()


# 5.1 Lie-group symmetry detection -----------------------------------------------------------------


reactions <- NULL
reactions <- addReaction(reactions, "TCA_buffer", "TCA_cell",
                         rate = "import_baso*TCA_buffer")
reactions <- addReaction(reactions, "TCA_cell", "TCA_buffer",
                         rate = "export_sinus*TCA_cell")
reactions <- addReaction(reactions, "TCA_cell", "TCA_cana",
                         rate = "export_cana*TCA_cell")
reactions <- addReaction(reactions, "TCA_cana", "TCA_buffer",
                         rate = "reflux*TCA_cana")
observables <- eqnvec(buffer = "s*TCA_buffer",
                      cellular = "s*(TCA_cana + TCA_cell)")
symmetryDetection(as.eqnvec(reactions), observables)

# 5.2 Analytical steady-state constraints ----------------------------------------------------------

steadyStates(reactions, file = "SS.Rds")


# Return back to old working directory

setwd(mywd)
