# Small reproducible community: 2 main levels x 2 nested levels x 5 samples
make_comm <- function(seed = 42) {
  set.seed(seed)
  main   <- rep(c("A", "B"), each = 10)
  nested <- rep(rep(c("x", "y"), each = 5), 2)
  lambda <- ifelse(nested == "x", 2, 6) + ifelse(main == "A", 0, 3)
  comm <- as.data.frame(sapply(1:8, function(j) stats::rpois(20, lambda * j / 4)))
  names(comm) <- paste0("sp", 1:8)
  list(comm = comm, main = main, nested = nested)
}
