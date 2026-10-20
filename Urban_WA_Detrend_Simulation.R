# ============================================================
# 只去掉叶面积长期趋势的条件 vine copula 对照
#
# 副主编要求的对照不是“每天都等于 1981–1990 均值”。
# 那个常数情景会同时去掉季节循环、年际波动和长期趋势。
#
# 本脚本的两个主情景：
#   Observed
#     观测到的城区、城郊叶面积。季节循环、年际波动和
#     长期趋势都保留。
#   Detrended
#     同一条序列只减去相对 1981–1990 的线性趋势。
#     季节循环和年际波动保留，序列在基线期的均值不变。
#
# 另外保留旧的 Constant 情景（每天都是 1981–1990 均值），
# 用来看原来的差值是不是去掉季节循环造成的。三个情景共用
# 同一个拟合好的 copula，只更换条件叶面积。
#
# 条件变量必须用“观测叶面积”的经验分布变到 (0,1)。
# 不能对常数序列或去趋势序列各自重新排秩，否则常数序列
# 会全部变成 0.5，去趋势带来的水平变化也会被排秩消掉。
#
# 不使用原来的全域谐波插值。那一步的预测只依赖年内日序，
# 趋势和年际异常在进入 copula 之前就已经没有了。
#
# 每个城市写出三份结果：
#   by_city/城市.csv          全时段平均和相关系数
#   daily/城市.rds           逐日 P、ET、P−ET，三情景各一列
#   annual/城市.csv           年尺度平均和年内相关系数
# 全部完成后合并为 city_summary.csv 和 annual_results.csv。
# ============================================================

args_all <- commandArgs(trailingOnly = TRUE)
n_cities_arg <- NA_integer_
n_cores_arg <- 1L
if (length(args_all) >= 1 && nzchar(args_all[1])) {
  n_cities_arg <- as.integer(args_all[1])
}
if (length(args_all) >= 2 && nzchar(args_all[2])) {
  n_cores_arg <- as.integer(args_all[2])
}

.libPaths(c(path.expand("~/Library/R/library"), .libPaths()))
suppressPackageStartupMessages({
  library(VineCopula)
  library(CDVineCopulaConditional)
  library(readxl)
})

BASE_DIR <- "/Volumes/202408/PminusE/AREA2018"
OUT_DIR <- file.path(BASE_DIR, "Vine_Detrend")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

START_DATE <- as.Date("1981-01-01")
END_DATE <- as.Date("2023-12-31")
dates <- seq.Date(START_DATE, END_DATE, by = "day")
stopifnot(length(dates) == 15705)
baseline <- format(dates, "%Y") %in% as.character(1981:1990)

# ------------------------------------------------------------
# 只去掉线性长期趋势。基线期中点处的趋势增量为 0，
# 所以 1981–1990 的平均水平保持不变。
# ------------------------------------------------------------
detrend_lai <- function(x, dates, baseline) {
  ok <- is.finite(x)
  out <- rep(NA_real_, length(x))
  if (sum(ok) < 365 || sum(ok & baseline) < 300) {
    return(x)
  }
  t <- as.numeric(dates)
  t0 <- mean(t[ok & baseline])
  fit <- lm(x[ok] ~ I(t[ok] - t0))
  slope <- unname(coef(fit)[2])
  out[ok] <- x[ok] - slope * (t[ok] - t0)
  attr(out, "slope_per_day") <- slope
  out
}

# 用参照样本（观测叶面积）的经验分布给任意新值一个伪观测。
# 与 VineCopula::pobs 的平均秩 / (n+1) 一致。
pobs_from_ref <- function(x_new, x_ref) {
  ref <- x_ref[is.finite(x_ref)]
  n <- length(ref)
  out <- rep(NA_real_, length(x_new))
  ok <- is.finite(x_new)
  if (n < 10 || !any(ok)) {
    return(out)
  }
  s <- sort(ref)
  le <- findInterval(x_new[ok], s)
  lt <- findInterval(x_new[ok], s, left.open = TRUE)
  u <- (lt + 0.5 * (le - lt)) / (n + 1)
  u <- pmin(pmax(u, 1 / (n + 1)), n / (n + 1))
  out[ok] <- u
  out
}

pobs_std <- function(x) {
  rank(x, na.last = "keep", ties.method = "average") / (sum(is.finite(x)) + 1)
}

pobs_zero <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  x_ok <- x[ok]
  zero_prop <- mean(x_ok == 0)
  if (zero_prop <= 0) {
    out[ok] <- rank(x_ok) / (length(x_ok) + 1)
    return(out)
  }
  idx_zero <- which(ok & x == 0)
  idx_pos <- which(ok & x > 0)
  if (length(idx_zero) > 0) {
    out[idx_zero] <- runif(length(idx_zero), 0, zero_prop * 0.99)
  }
  if (length(idx_pos) > 0) {
    r <- rank(x[idx_pos]) / (length(idx_pos) + 1)
    out[idx_pos] <- zero_prop + (1 - zero_prop) * r
  }
  out
}

inv_precip <- function(u, p_hist) {
  u[is.na(u)] <- 0.5
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  x_pos <- p_hist[is.finite(p_hist) & p_hist > 0]
  if (length(x_pos) == 0) {
    return(rep(0, length(u)))
  }
  zero_prop <- mean(p_hist == 0, na.rm = TRUE)
  out <- rep(0, length(u))
  idx <- which(u > zero_prop)
  if (length(idx) > 0) {
    uu <- (u[idx] - zero_prop) / (1 - zero_prop)
    uu <- pmin(pmax(uu, 1e-6), 1 - 1e-6)
    out[idx] <- as.numeric(quantile(x_pos, uu, type = 8, names = FALSE))
  }
  out
}

inv_std <- function(u, x_hist) {
  u[is.na(u)] <- 0.5
  u <- pmin(pmax(u, 1e-6), 1 - 1e-6)
  x_ok <- x_hist[is.finite(x_hist)]
  as.numeric(quantile(x_ok, u, type = 8, names = FALSE))
}

read_city_matrix <- function(path) {
  raw <- read_excel(path, col_names = TRUE, progress = FALSE)
  date_col <- raw[[1]]
  mat <- as.matrix(raw[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  if (nrow(mat) != length(dates)) {
    stop(sprintf("%s 有 %d 行，期望 %d 行。", path, nrow(mat), length(dates)))
  }
  list(cities = colnames(mat), values = mat, date_head = as.character(date_col[1]))
}

load_inputs <- function() {
  cache <- file.path(OUT_DIR, "inputs_cache.rds")
  if (file.exists(cache)) {
    message("读取缓存 ", cache)
    return(readRDS(cache))
  }
  message("首次读取 Excel，完成后写入缓存。")
  pieces <- list(
    P = read_city_matrix(file.path(BASE_DIR, "P19812023UR.xlsx")),
    ET = read_city_matrix(file.path(BASE_DIR, "ET19812023UR.xlsx")),
    LAIu = read_city_matrix(file.path(BASE_DIR, "LAI19812023ur.xlsx")),
    LAIp = read_city_matrix(file.path(BASE_DIR, "LAI19812023vi.xlsx"))
  )
  cities <- pieces$P$cities
  for (nm in names(pieces)) {
    if (!setequal(pieces[[nm]]$cities, cities)) {
      stop(nm, " 的城市列与降水不一致。")
    }
  }
  order_to <- function(obj) obj$values[, cities, drop = FALSE]
  out <- list(
    cities = cities,
    P = order_to(pieces$P),
    ET = order_to(pieces$ET),
    LAIu = order_to(pieces$LAIu),
    LAIp = order_to(pieces$LAIp)
  )
  saveRDS(out, cache, compress = "xz")
  out
}

# D-vine 的条件列必须按 RVM$Matrix 对角线的变量编号排列。
condition_matrix <- function(U, RVM, nx = 2) {
  idx <- diag(RVM$Matrix)[seq_len(nx)]
  U[, idx, drop = FALSE]
}

simulate_one <- function(city_i, inputs) {
  city <- inputs$cities[city_i]
  seed_case <- 1000L + city_i

  P <- inputs$P[, city_i]
  ET <- inputs$ET[, city_i]
  LAIu <- inputs$LAIu[, city_i]
  LAIp <- inputs$LAIp[, city_i]

  LAIu_dt <- detrend_lai(LAIu, dates, baseline)
  LAIp_dt <- detrend_lai(LAIp, dates, baseline)
  LAIu_const <- rep(mean(LAIu[baseline], na.rm = TRUE), length(dates))
  LAIp_const <- rep(mean(LAIp[baseline], na.rm = TRUE), length(dates))

  # 零降水的随机化只服务边缘分布，不进入三个情景的比较。
  set.seed(seed_case)
  u_P <- pobs_zero(P)
  u_ET <- pobs_std(ET)
  u_Lu <- pobs_std(LAIu)
  u_Lp <- pobs_std(LAIp)

  # 去趋势值和常数都映射到观测叶面积的分布上。
  u_Lu_dt <- pobs_from_ref(LAIu_dt, LAIu)
  u_Lp_dt <- pobs_from_ref(LAIp_dt, LAIp)
  u_Lu_c <- pobs_from_ref(LAIu_const, LAIu)
  u_Lp_c <- pobs_from_ref(LAIp_const, LAIp)

  U <- cbind(u_P, u_ET, u_Lu, u_Lp)
  ok <- stats::complete.cases(U)
  if (sum(ok) < 300) {
    return(NULL)
  }

  RVM <- CDVineCondFit(
    U[ok, , drop = FALSE],
    Nx = 2,
    type = "DVine",
    treecrit = "BIC",
    selectioncrit = "BIC"
  )

  run_case <- function(u1, u2) {
    UU <- cbind(u_P, u_ET, u1, u2)
    cond <- condition_matrix(UU, RVM, nx = 2)
    sim <- CDVineCondSim(RVM, cond)
    # 返回矩阵前两列是被模拟的变量，后两列是条件变量。
    P_sim <- inv_precip(sim[, 1], P)
    ET_sim <- inv_std(sim[, 2], ET)
    list(P = P_sim, ET = ET_sim, WA = P_sim - ET_sim)
  }

  # 三个情景重置到同一个种子，潜在随机数相同，只更换叶面积条件。
  set.seed(seed_case)
  obs <- run_case(u_Lu, u_Lp)
  set.seed(seed_case)
  det <- run_case(u_Lu_dt, u_Lp_dt)
  set.seed(seed_case)
  con <- run_case(u_Lu_c, u_Lp_c)

  pack <- function(s) {
    c(
      P = mean(s$P, na.rm = TRUE),
      ET = mean(s$ET, na.rm = TRUE),
      WA = mean(s$WA, na.rm = TRUE),
      corr = suppressWarnings(cor(s$P, s$ET, use = "complete.obs"))
    )
  }
  o <- pack(obs)
  d <- pack(det)
  cnst <- pack(con)

  daily <- cbind(
    P_obs = obs$P, ET_obs = obs$ET, WA_obs = obs$WA,
    P_det = det$P, ET_det = det$ET, WA_det = det$WA,
    P_const = con$P, ET_const = con$ET, WA_const = con$WA
  )
  # 行顺序与 dates 一致：1981-01-01 到 2023-12-31，共 15705 天。
  rownames(daily) <- as.character(dates)

  summary_row <- data.frame(
    city = city,
    n_days = sum(ok),
    lai_u_slope_per_year = attr(LAIu_dt, "slope_per_day") * 365.25,
    lai_p_slope_per_year = attr(LAIp_dt, "slope_per_day") * 365.25,
    P_obs = o["P"], ET_obs = o["ET"], WA_obs = o["WA"], corr_obs = o["corr"],
    P_det = d["P"], ET_det = d["ET"], WA_det = d["WA"], corr_det = d["corr"],
    P_const = cnst["P"], ET_const = cnst["ET"], WA_const = cnst["WA"], corr_const = cnst["corr"],
    dWA_trend = unname(o["WA"] - d["WA"]),
    dCorr_trend = unname(o["corr"] - d["corr"]),
    dP_trend = unname(o["P"] - d["P"]),
    dET_trend = unname(o["ET"] - d["ET"]),
    dWA_const = unname(o["WA"] - cnst["WA"]),
    dCorr_const = unname(o["corr"] - cnst["corr"]),
    stringsAsFactors = FALSE
  )
  list(summary = summary_row, daily = daily, annual = annual_from_daily(city, daily))
}

annual_from_daily <- function(city, daily) {
  yr <- as.integer(format(dates, "%Y"))
  rows <- lapply(sort(unique(yr)), function(y) {
    b <- daily[yr == y, , drop = FALSE]
    m <- colMeans(b, na.rm = TRUE)
    data.frame(
      city = city,
      year = y,
      P_obs = m[["P_obs"]], ET_obs = m[["ET_obs"]], WA_obs = m[["WA_obs"]],
      P_det = m[["P_det"]], ET_det = m[["ET_det"]], WA_det = m[["WA_det"]],
      P_const = m[["P_const"]], ET_const = m[["ET_const"]], WA_const = m[["WA_const"]],
      corr_obs = suppressWarnings(cor(b[, "P_obs"], b[, "ET_obs"], use = "complete.obs")),
      corr_det = suppressWarnings(cor(b[, "P_det"], b[, "ET_det"], use = "complete.obs")),
      corr_const = suppressWarnings(cor(b[, "P_const"], b[, "ET_const"], use = "complete.obs")),
      dWA_trend = unname(m[["WA_obs"]] - m[["WA_det"]]),
      dWA_const = unname(m[["WA_obs"]] - m[["WA_const"]]),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

out_subdir <- function(name) {
  d <- file.path(OUT_DIR, name)
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  d
}

save_city_outputs <- function(result) {
  city <- result$summary$city
  write.csv(result$summary, file.path(out_subdir("by_city"), paste0(city, ".csv")), row.names = FALSE)
  saveRDS(result$daily, file.path(out_subdir("daily"), paste0(city, ".rds")), compress = FALSE)
  write.csv(result$annual, file.path(out_subdir("annual"), paste0(city, ".csv")), row.names = FALSE)
}

done_cities <- function() {
  daily_files <- list.files(file.path(OUT_DIR, "daily"), pattern = "\\.rds$", full.names = FALSE)
  annual_files <- list.files(file.path(OUT_DIR, "annual"), pattern = "\\.csv$", full.names = FALSE)
  summary_files <- list.files(file.path(OUT_DIR, "by_city"), pattern = "\\.csv$", full.names = FALSE)
  if (!length(daily_files) || !length(annual_files) || !length(summary_files)) {
    return(character())
  }
  daily_ids <- tools::file_path_sans_ext(daily_files)
  annual_ids <- tools::file_path_sans_ext(annual_files)
  summary_ids <- tools::file_path_sans_ext(summary_files)
  Reduce(intersect, list(daily_ids, annual_ids, summary_ids))
}

combine_summaries <- function() {
  files <- list.files(file.path(OUT_DIR, "by_city"), pattern = "\\.csv$", full.names = TRUE)
  if (!length(files)) {
    return(invisible(NULL))
  }
  out <- do.call(rbind, lapply(files, utils::read.csv, stringsAsFactors = FALSE))
  out <- out[!duplicated(out$city), , drop = FALSE]
  utils::write.csv(out, file.path(OUT_DIR, "city_summary.csv"), row.names = FALSE)
  annual_files <- list.files(file.path(OUT_DIR, "annual"), pattern = "\\.csv$", full.names = TRUE)
  if (length(annual_files)) {
    annual <- do.call(rbind, lapply(annual_files, utils::read.csv, stringsAsFactors = FALSE))
    utils::write.csv(annual, file.path(OUT_DIR, "annual_results.csv"), row.names = FALSE)
  }
  out
}

main <- function() {
  inputs <- load_inputs()
  cities <- inputs$cities
  if (!is.na(n_cities_arg) && n_cities_arg > 0) {
    cities <- cities[seq_len(min(n_cities_arg, length(cities)))]
  }
  done <- done_cities()
  todo <- which(inputs$cities %in% setdiff(cities, done))
  message(sprintf("待计算 %d 个城市，已完成 %d 个。", length(todo), length(done)))
  if (length(todo) == 0) {
    return(invisible(NULL))
  }

  one <- function(i) {
    t0 <- proc.time()[3]
    row <- tryCatch(
      simulate_one(i, inputs),
      error = function(e) {
        message("失败 ", inputs$cities[i], ": ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(row)) {
      save_city_outputs(row)
      dt <- proc.time()[3] - t0
      message(sprintf(
        "%s  dWA_trend=%.4f  dWA_const=%.4f  用时 %.1f 秒",
        row$summary$city, row$summary$dWA_trend, row$summary$dWA_const, dt
      ))
    }
    row$summary
  }

  if (n_cores_arg > 1 && length(todo) > 1) {
    parallel::mclapply(todo, one, mc.cores = n_cores_arg, mc.preschedule = FALSE)
  } else {
    for (i in todo) one(i)
  }
  tab <- combine_summaries()
  if (!is.null(tab)) {
    message(sprintf(
      "完成 %d 个城市。去趋势对照：ΔWA 均值 %.4f（正值比例 %.1f%%）；常数对照：ΔWA 均值 %.4f（正值比例 %.1f%%）。",
      nrow(tab),
      mean(tab$dWA_trend, na.rm = TRUE),
      100 * mean(tab$dWA_trend > 0, na.rm = TRUE),
      mean(tab$dWA_const, na.rm = TRUE),
      100 * mean(tab$dWA_const > 0, na.rm = TRUE)
    ))
  }
  message("写出 ", file.path(OUT_DIR, "city_summary.csv"))
}

main()
