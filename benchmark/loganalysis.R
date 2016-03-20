library(tcR)

em <- list()
em.data <- list()
sg <- list()
sg.data <- list()
setwd("~/Projects/ymir/benchmark/")
sg.df <- matrix(0, 0, 5)
for (file in list.files("./log", full.names = T)) {
  spl <- strsplit(file, "_", T)[[1]]
  # vec <- list(as.numeric(readLines(file))[1:9])
  vec <- list(as.numeric(readLines(file)))
  if (tail(spl, 1) == "1.ll.txt") {
    if (spl[2] == "em") {
      em <- c(em, list(as.numeric(c(spl[4]))))
      em.data <- c(em.data, vec)
    } else {
      sg <- c(sg, list(as.numeric(c(spl[4], spl[6], spl[8], spl[10], spl[12]))))
      sg.data <- c(sg.data, vec)
      sg.df <- rbind(sg.df, c(vec[[1]][9], as.numeric(c(spl[6], spl[8], spl[10], spl[12]))))
    }
  }
}

sg.df <- as.data.frame.matrix(sg.df)
colnames(sg.df) <- c("logL", "Block", "Alpha", "Beta", "K")
maxB1 = filter(group_by(sg.df, Block), logL == max(logL))
maxB2 = filter(group_by(sg.df, Block), logL < max(logL))
maxB2 = filter(maxB2, logL == max(logL))

remove.iters <- c(1:1)

sg <- as.data.frame(t(as.data.frame(sg)))
sg.data <- as.data.frame(sg.data)[-remove.iters, ]
colnames(sg) <- c("niter", "block", "alpha", "beta", "K")
colnames(sg.data) <- paste0("block_", sg$block, "_alpha_", sg$alpha, "_beta_", sg$beta, "_K_", sg$K)
sg.data$Iteration <- 1:nrow(sg.data)
sg.data$Processed <- sg.data$block
sg.data <- melt(sg.data, id.vars = "Iteration")
names(sg.data)[2:3] <- c("Algorithm", "logL")
sg.data$Data <- sg.data$Iteration * as.numeric(sapply(strsplit(as.character(sg.data$Algorithm), "_", T), "[[", 2))
sg.data$Block <- sapply(strsplit(as.character(sg.data$Algorithm), "_", T), "[[", 2)

em.data <- data.frame(logL = em.data[[1]][-remove.iters])
em.data$Iteration <- 1:nrow(em.data)
em.data$Algorithm <- "EM"
em.data$Data <- em.data$Iteration * 100000

p1 <- ggplot() +
  geom_line(aes(x = Iteration, y = logL, colour = Algorithm, group = Algorithm), linetype = "dashed", size = .5, data = sg.data) +
  geom_point(aes(x = Iteration, y = logL, colour = Algorithm, fill = Algorithm, shape = Block), size = 4, data = sg.data) +
  geom_line(aes(x = Iteration, y = logL, colour = Algorithm, group = Algorithm), size = .5, data = em.data) +
  geom_point(aes(x = Iteration, y = logL, colour = Algorithm, fill = Algorithm, shape = "EM"), size = 4, data = em.data) +
  xlab("Итерация") + ggtitle("A") +
  scale_colour_discrete(guide = guide_legend(title = "Алгоритм")) + 
  scale_fill_discrete(guide = guide_legend(title = "Алгоритм")) + 
  scale_shape_discrete(guide = guide_legend(title = "m")) + 
  theme_linedraw()

p2 <- ggplot() +
  geom_line(aes(x = Data, y = logL, colour = Algorithm, group = Algorithm), linetype = "dashed", size = .5, data = sg.data) +
  geom_point(aes(x = Data, y = logL, colour = Algorithm, fill = Algorithm, shape = Block), size = 4, data = sg.data) +
  geom_line(aes(x = Data, y = logL, colour = Algorithm, group = Algorithm), size = .5, data = em.data) +
  geom_point(aes(x = Data, y = logL, colour = Algorithm, fill = Algorithm, shape = "EM"), size = 4, data = em.data) +
  scale_x_log10() +
  xlab("Обработанные объекты") + ggtitle("B") + 
  scale_colour_discrete(guide = guide_legend(title = "Алгоритм")) + 
  scale_fill_discrete(guide = guide_legend(title = "Алгоритм")) + 
  scale_shape_discrete(guide = guide_legend(title = "m")) + 
  theme_linedraw()

grid.arrange(p1, p2, ncol = 1)

# em.pars <- c()
# sg.pars <- list()
# for (file in list.files("./log", full.names = T)) {
#   spl <- strsplit(file, "_", T)[[1]]
#   # vec <- list(as.numeric(readLines(file))[1:9])
#   vec <- as.numeric(readLines(file))
#   if (tail(spl, 1) == "1.pars.txt") {
#     if (spl[2] == "em") {
#       em.pars <- vec
#     } else {
#       sg.pars <- c(sg.pars, list(vec))
#     }
#   }
# }
# 
# cor.names <- paste0("m=", sg$block, ";a=", sg$alpha, ";b=", sg$beta, ";k=", sg$K)
# cor.plots <- lapply(1:length(sg.pars), function (j) { 
#   df <- data.frame(EM = em.pars, 
#                    SG = sg.pars[[j]])
#   ggplot() + geom_point(aes(x = EM, y = SG), data = df) + 
#     geom_text(aes(x = .85, y = .15, size = 1, colour = "red", label = round(cor(em.pars, sg.pars[[j]]), 2))) +
#     scale_colour_discrete(guide=FALSE) +
#     scale_size_continuous(guide=FALSE) +
#     theme_linedraw() +
#     geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#     coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
#     ggtitle(cor.names[j]) + xlab("") + ylab("")
#   })
# do.call(grid.arrange, c(cor.plots, list(ncol = 5)))