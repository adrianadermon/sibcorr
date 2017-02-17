# This file simulates a sibling dataset and tests some code for estimating the
# sibling correlation using the method described in Solon, Page, and Duncan (2000)
# and applied by Hällsten (2014) on Swedish data.
# =====================================================================

library(feather)
library(data.table)
library(sibcorr)

df <- read_feather("../../data_small.feather")

# Rename columns
names(df) <- c("id_mormor", "id_mor", "id_barn", "gender", "by", "w")

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

context("Correlation estimation")

test_that("sibling correlations are estimated correctly", {
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df),
    c(correlation = 0.47941563781226687, n_individuals = 39282, n_pairs = 39352, n_families = 11645))
  expect_equal(
    sibcorr(w ~ factor(by) + factor(gender) | id_barn + id_mor, data = df),
    c(correlation = 0.49059087767216614, n_individuals = 39282, n_pairs = 39352, n_families = 11645))
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, weight = 1),
    c(correlation = 0.47876571240015514, n_individuals = 39282, n_pairs = 39352, n_families = 11645))
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, weight = 2),
    c(correlation = 0.47966007916912057, n_individuals = 39282, n_pairs = 39352, n_families = 11645))
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, weight = 3),
    c(correlation = 0.47991978060871299, n_individuals = 39282, n_pairs = 39352, n_families = 11645))
})

test_that("cousin correlations are estimated correctly", {
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor + id_mormor, data = df),
    c(correlation = 0.139953, n_individuals = 39282, n_pairs = 76208, n_families = 5017),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ factor(by) + factor(gender) | id_barn + id_mor + id_mormor, data = df),
    c(correlation = 0.1440984, n_individuals = 39282, n_pairs = 76208, n_families = 5017),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor + id_mormor, data = df, weight = 1),
    c(correlation = 0.1376879, n_individuals = 39282, n_pairs = 76208, n_families = 5017),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor + id_mormor, data = df, weight = 2),
    c(correlation = 0.1386408, n_individuals = 39282, n_pairs = 76208, n_families = 5017),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor + id_mormor, data = df, weight = 3),
    c(correlation = 0.1387906, n_individuals = 39282, n_pairs = 76208, n_families = 5017),
    tolerance = 0.00001)
})


test_that("restrictions work correctly", {
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, restriction = c("by", 2)),
    c(correlation = 0.4719575, n_individuals = 39282, n_pairs = 3511, n_families = 2750),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, restriction = c("by", "unequal")),
    c(correlation = 0.4769372, n_individuals = 39282, n_pairs = 37466, n_families = 11380),
    tolerance = 0.00001)
  expect_equal(
    sibcorr(w ~ 0 | id_barn + id_mor, data = df, restriction = c("by", 3, 6)),
    c(correlation = 0.4821474, n_individuals = 39282, n_pairs = 12177, n_families = 6330),
    tolerance = 0.00001)
  expect_equal(
    unname(sibcorr(w ~ factor(by) + factor(gender) | id_barn + id_mor, data = df, restriction = c("by", 5, 2))),
    c(0.490645554813508, 39282, 26500, 9786))
  expect_equal(
    unname(sibcorr(w ~ 0 | id_barn + id_mor + id_mormor, data = df, restriction = c("by", 0, 9))),
    c(0.143199878287237, 39282, 55119, 4871))
  expect_equal(
    unname(sibcorr(w ~ factor(by) + factor(gender) | id_barn + id_mor + id_mormor, data = df, restriction = c("gender", "unequal"))),
    c(0.145962732074221, 39282, 38069, 4573))
})

context("Bootstrap estimation")

set.seed(20170206)

test_that("bootstrap works correctly", {
  expect_equal(
    sibcorr_bs(w ~ 0 | id_barn + id_mor, data = df),
    c(correlation = 0.47941563781226687, n_individuals = 39282, n_pairs = 39352, n_families = 11645, "2.5%" = 0.46699407161779294, "97.5%" = 0.49040559449898979))
  expect_equal(
    sibcorr_bs(w ~ 0 | id_barn + id_mor + id_mormor, data = df),
    c(correlation = 0.1399530, n_individuals = 39282, n_pairs = 76208, n_families = 5017, "2.5%" = 0.1278036, "97.5%" = 0.1494802),
    tolerance = 0.0000001)
})

set.seed(20170211)

test_that("multi-processor bootstrap works correctly", {
  expect_equal(
    sibcorr_bs(w ~ 0 | id_barn + id_mor, data = df, cores = 4),
    c(correlation = 0.47941563781226687, n_individuals = 39282, n_pairs = 39352, n_families = 11645, "2.5%" = 0.47, "97.5%" = 0.49),
    tolerance = 0.01)
  expect_equal(
    sibcorr_bs(w ~ 0 | id_barn + id_mor + id_mormor, data = df, cores = 4),
    c(correlation = 0.13978219591128127, n_individuals = 39282, n_pairs = 76208, n_families = 5017, "2.5%" = 0.13, "97.5%" = 0.15),
    tolerance = 0.01)
})

set.seed(20170217)

test_that("bootstrap works with restriction", {
  expect_equal(
    unname(sibcorr_bs(w ~ 0 | id_barn + id_mor, data = df, restriction = c("by", 4))),
    c(0.499471797184667, 39282, 3029, 2440, 0.465584526329792, 0.52979876204937))
  expect_equal(
    unname(sibcorr_bs(w ~ 0 | id_barn + id_mor, data = df, restriction = c("by", 4), cores = 4)),
    c(0.499471797184667, 39282, 3029, 2440, 0.47, 0.54),
    tolerance = 0.01)
})
