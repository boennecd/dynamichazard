#####
# Load data example to play arround with and validate against
is_censored = c(6, 27, 34, 36, 42, 46, 48:51,
                51 + c(15, 30:28, 33, 35:37, 39, 40, 42:45))
head_neck_cancer = data.frame(
  id = 1:96,
  start = rep(0, 96),
  stop = c(
    1, 2, 2, rep(3, 6), 4, 4, rep(5, 8),
    rep(6, 7), 7, 8, 8, 8,
    9, 9, 10, 10, 10, 11, 14, 14, 14, 15, 18, 18, 20, 20, 37,
    37, 38, 41, 45, 47, 47,

    2, 2, 3, rep(4, 4), rep(5, 5), rep(6, 5),
    7, 7, 7, 9, 10, 11, 12, 15, 16, 18, 18, 18, 21,
    21, 24, 25, 27, 36, 41, 44, 52, 54, 59, 59, 63, 67, 71, 76),
  event = !(1:96 %in% is_censored),
  group = factor(c(rep(1, 45 + 6), rep(2, 45))))

rm(is_censored)

head_neck_cancer$group = factor(head_neck_cancer$group, levels = c(2, 1))
