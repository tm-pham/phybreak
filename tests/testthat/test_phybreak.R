library(phybreak)
context("phybreak construction")

test_that("use.tree works correctly", {
  set.seed(0)
  sim1 <- sim_phybreak(5)
  pbo1 <- phybreak(sim1, use.tree = TRUE)

  expect_equal(pbo1$v$infectors,
               match(sim1$sim.infectors, c("index", paste0("host.", 1:5))) - 1)
  expect_equal(pbo1$d$reference.date + pbo1$v$inftimes,
               unname(sim1$sim.infection.times[pbo1$d$hostnames]))
  expect_equal(pbo1$d$reference.date, 
               min(sim1$sample.times))
  expect_equal(pbo1$d$reference.date + pbo1$v$nodetimes[1:9], 
               sim1$sim.tree$root.edge + ape::node.depth.edgelength(sim1$sim.tree))
})
