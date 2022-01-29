ind <- matrix(rnorm(670000*132), 670000, 132)
ksvd <- svd(cor(ind))
x <- cor(ind)
my_e <- eigen(x,only.values = TRUE, symmetric = TRUE)
? cor

my_e <- eigen(t(ind) %*% (ind),only.values = TRUE, symmetric = FALSE)

ind_z <- scale(ind)
my_e <- eigen(t(ind_z) %*% (ind_z),only.values = TRUE, symmetric = TRUE)
sessionInfo()

a <- Sys.time()
x <- cor(ind)
Sys.time()-a # 7.23

a <- Sys.time()
x2 <- crossprod(scale(ind, TRUE, TRUE))/(nrow(x)-1)
Sys.time()-a #10.26

install.packages("coop")
library(coop)
a <- Sys.time()
x2 <- coop::pcor(ind)
Sys.time()-a #4.73
max(abs(x - x2))

sessionInfo()
# install openblas
sudo apt-get install libopenblas-base
sudo update-alternatives --config libblas.so.3-aarch64-linux-gnu

kevin <- 3.14

