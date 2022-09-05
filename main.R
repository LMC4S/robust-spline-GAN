# Running jobs using shell scripts (Rscript), call main.R in terminal with
# ...arguments.

# can also run directly in R with default parameters

rm(list=ls())
options(warn = 1)

library(tictoc)
suppressMessages(library(CVXR))
suppressMessages(library(tidyverse))

setwd("~/robust-spline-GAN")
source('source/utilities.R')
source('source/fGAN.R')
source('source/train.R')

getwd()

cmd_args <- commandArgs(trailingOnly=TRUE)

defaults <- list( 
  eps = 0.1, # contamination proportion
  p = 5, # dimension
  n = 1000, # sample size
  iter = 50, # iteration times
  repeat_times = 1, # fold of repeated experiments
  Loss = "rKL", # Loss, "JS", "rKL", or "Hinge"
  optim = "MOSEK", # convex solver. "MOSEK" or "ECOS".
  lr = 1, # generator learning rate
  pen = 0.3, # penalty level. Use 0.025 for JS, 0.3 for rKL, 0.1 for hinge 
  Q = 'A', # Type of contamination, "A" or "B".
  l1 = TRUE, # l1 penalty or l2, logical.
  wd= "./results", # output directory
  seed = 0 # random seed
  ) ## default values of training and model setting arguments 

if (length(defaults) == length(cmd_args)) {
  for (i in 1:length(defaults)) {
    if (class(defaults[[i]]) == "numeric") {
      defaults[[i]]  <- as.numeric(cmd_args[i])
    } else if (cmd_args[[i]] %in% c("FALSE", "TRUE", "T", "F")) {
      defaults[[i]]  <- as.logical(cmd_args[i])
    } else {
      defaults[[i]] <- cmd_args[i]
    }
    print(defaults[[i]])
  }
} else {
  warning("Argument missing. Using default for all arguments")
}
list2env(defaults, globalenv())

## parse command args, load into global environment
## if any variable named in defaults doesn't exist, then create it
## with value from defaults


args <- list(eps=eps, p=p, n=n, Cov="Ar", sd=1, mu=NULL, # model spec
             Q=Q, # contamination setting, Q = "A" or "B"
             fold=repeat_times,
             pen_lvl=pen,
             L1_pen=l1,
             batch_x=n, batch_z=n,
             kendall_size=n,
             lr=lr,
             decay_lr=T,
             iter=iter,
             optim='MOSEK', # MOSEK or ECOS solver
             Loss=Loss, # "rKL", "JS", or "Hinge"
             seed=seed)


print(paste("\n\n Task seed is:", args$seed))

dir.create(file.path(wd))
setwd(file.path(wd))

sink(paste("seed", args$seed,  "_",  Loss, "_p", p, "_n", n, "_eps",
           100*eps, "_pen", pen, "_l1", l1,
          "_Qtype", Q,
           ".log", sep=""))

Sys.time()

tic()
res <- do.call(train, args=args)
cat("\n Running time: ")
toc()

rm(list=setdiff(ls(), c("res", "args", "cmd_args")))
save.image(paste("seed", args$seed, "_", args$Loss, "_p", args$p, "_n", args$n, "_eps",
                 100*args$eps, "_pen", args$pen_lvl, "_l1", args$L1_pen,
                 "_Qtype", args$Q,
                 ".RData", sep=""))

if (args$fold != 1) { #print out result summary if it's repeated simulation
    summary(res$L2)
    summary(res$fro)
    summary(res$op)
}

print(args)
sessionInfo()

