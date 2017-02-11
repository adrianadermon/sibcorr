# This file simulates a sibling dataset and tests some code for estimating the
# sibling correlation using the method described in Solon, Page, and Duncan (2000)
# and applied by Hällsten (2014) on Swedish data.
# =====================================================================

library(feather)
library(data.table)
library(sibcorr)

context("Correlation estimation")


df <- read_feather("../../data_small.feather")

# Convert columns to factors
df$gender <- as.factor(df$gender)
df$by <- as.factor(df$by)

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

test_that("sibling correlations are estimated correctly", {
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor"),
    0.47941563781226687)
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", controls = c("by", "gender")),
    0.49059087767216614)
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", weight = 1),
    0.47876571240015514)
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", weight = 2),
    0.47966007916912057)
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", weight = 3),
    0.47991978060871299)
})

test_that("cousin correlations are estimated correctly", {
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", cousins = TRUE),
    0.13978219591128127)
  expect_equal(
    sibcorr(df, y = "w", id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", controls = c("by", "gender"), cousins = TRUE),
    0.14397525365735642)
  expect_equal(
    sibcorr(df, y = "w", "id_barn", "id_mor", "id_mormor", weight = 1, cousins = TRUE),
    0.13751981943694128)
  expect_equal(
    sibcorr(df, y = "w", "id_barn", "id_mor", "id_mormor", weight = 2, cousins = TRUE),
    0.13847156469002828)
  expect_equal(
    sibcorr(df, y = "w", "id_barn", "id_mor", "id_mormor", weight = 3, cousins = TRUE),
    0.13862117315575975)
})

set.seed(20170206)

test_that("bootstrap works correctly", {
  expect_equal(
    sibcorr_bs(df, y = "w", id1 = "id_barn", id2 = "id_mor"),
    c(0.47941563781226687, "2.5%" = 0.46699407161779294, "97.5%" = 0.49040559449898979))
  expect_equal(
    sibcorr_bs(df, y = "w", id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", cousins = TRUE),
    c(0.13978219591128127, "2.5%" = 0.12798557949883563, "97.5%" = 0.14930428951834102))
})

set.seed(20170211)

test_that("multi-processor bootstrap works correctly", {
  expect_equal(
    sibcorr_bs(df, y = "w", id1 = "id_barn", id2 = "id_mor", cores = 4),
    c(0.47941563781226687, "2.5%" = 0.47, "97.5%" = 0.49),
    tolerance = 0.01)
  expect_equal(
    sibcorr_bs(df, y = "w", id1 = "id_barn", id2 = "id_mor", id3 = "id_mormor", cousins = TRUE, cores = 4),
    c(0.13978219591128127, "2.5%" = 0.13, "97.5%" = 0.15),
    tolerance = 0.01)
})
