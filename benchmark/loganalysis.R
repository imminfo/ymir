library(ggplot2)

em <- list()
em.data <- list()
sg <- list()
sg.data <- list()
setwd("~/Projects/ymir/benchmark/")
for (file in list.files("./log", full.names = T)) {
  spl <- strsplit(file, "_", T)[[1]]
  vec <- list(as.numeric(readLines(file)))
  if (spl[2] == "em") {
    em <- c(em, list(as.numeric(c(spl[4]))))
    em.data <- c(em.data, vec)
  } else {
    sg <- c(sg, list(as.numeric(c(spl[4], spl[6], spl[8], spl[10]))))
    sg.data <- c(sg.data, vec)
  }
}

remove.iters <- c(1:1)

sg <- as.data.frame(t(as.data.frame(sg)))
sg.data <- as.data.frame(sg.data)[-remove.iters, ]
colnames(sg) <- c("niter", "block", "alpha", "K")
colnames(sg.data) <- paste0("block_", sg$block, "_alpha_", sg$alpha, "_K_", sg$K)
sg.data$Iteration <- 1:nrow(sg.data)
sg.data$Processed <- sg.data$block 
sg.data <- melt(sg.data, id.vars = "Iteration")
names(sg.data)[2:3] <- c("Algorithm", "logL")
sg.data$Data <- sg.data$Iteration * as.numeric(sapply(strsplit(as.character(sg.data$Algorithm), "_", T), "[[", 2))


em.data <- data.frame(logL = em.data[[1]][-remove.iters])
em.data$Iteration <- 1:nrow(em.data)
em.data$Algorithm <- "EM"
em.data$Data <- em.data$Iteration * 10000

p1 <- ggplot() +
  geom_line(aes(x = Iteration, y = logL, colour = Algorithm, group = Algorithm), size = .8, linetype = 2, data = sg.data) +
  geom_line(aes(x = Iteration, y = logL, colour = Algorithm, group = Algorithm), size = .8, data = em.data) +
  theme_linedraw()

p2 <- ggplot() +
  geom_line(aes(x = Data, y = logL, colour = Algorithm, group = Algorithm), size = .8, linetype = 2, data = sg.data) +
  geom_line(aes(x = Data, y = logL, colour = Algorithm, group = Algorithm), size = .8, data = em.data) +
  scale_x_log10() +
  theme_linedraw()

grid.arrange(p1, p2, ncol = 1)