# This file simulates a sibling dataset and tests some code for estimating the
# sibling correlation using the method described in Solon, Page, and Duncan (2000)
# and applied by Hällsten (2014) on Swedish data.
# =====================================================================

library(feather)
library(data.table)
library(sibcorr)

context("Correlation estimation")


dt <- read_feather("../../data_small.feather")

# Convert to data table and set key
setDT(dt, key = "id2")

# Convert columns to factors
for (col in c("by", "gender")) dt[ , (col) := as.factor(dt[[col]])]

setnames(dt, "id1", "id_barn")
setnames(dt, "id2", "id_mor")
setnames(dt, "id3", "id_mormor")

### Required naming scheme:
# Outcome variable = y
# Individual identifier = id1
# Family identifier = id2
# Extended family (cousin group) identifier = id3

### Syntax:
# sibcorr(data, weight = 4, controls = NULL, cousins = FALSE)
# data should be a data frame or data table
# weight selects the weighting scheme among these options (default is 4):
# 1 (all families weighted equally, strongest down-weighting of large families) = [1/2n(n - 1)]^{-1}
# 2 (strong down-weighting) = [1/2(n - 1)]^{-1}
# 3 (weak down-weighting) = [1/2[n(n - 1)]^{1/2}]^{-1}
# 4 (equal cousin pair weight, no down-weighting of large families) = 1
# controls should be a vector of control variables that are to be regressed out before estimating the correlation
# cousins selects estimation of the sibling or cousin correlation. Set cousins = TRUE for cousin correlation

test_that("sibling correlations are estimated correctly", {
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor"), 0.47941563781226687)
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", controls = c("by", "gender")), 0.49059087767216614)
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", weight = 1), 0.47876571240015514)
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", weight = 2), 0.47966007916912057)
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", weight = 3), 0.47991978060871299)
})

test_that("cousin correlations are estimated correctly", {
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", cousins = TRUE), 0.13978219591128127)
  expect_equal(sibcorr(dt, id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", controls = c("by", "gender"), cousins = TRUE), 0.14397525365735642)
  expect_equal(sibcorr(dt, "id_barn", "id_mor", "id_mormor", weight = 1, cousins = TRUE), 0.13751981943694128)
  expect_equal(sibcorr(dt, "id_barn", "id_mor", "id_mormor", weight = 2, cousins = TRUE), 0.13847156469002828)
  expect_equal(sibcorr(dt, "id_barn", "id_mor", "id_mormor", weight = 3, cousins = TRUE), 0.13862117315575975)
})



# Estimate sibling correlation with boostrap standard errors
#sibcorr_bs(dt, controls = c("by", "gender"), cousins = FALSE, reps = 10, weight = 2)

# Estimate cousin correlation with boostrap standard errors
#sibcorr_bs(dt, controls = c("by", "gender"), cousins = TRUE, reps = 10, weight = 2)
