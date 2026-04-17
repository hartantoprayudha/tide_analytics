#!/usr/bin/env Rscript

# =============================================================
# Script Analisis Pasang Surut (Pasut)
# Fitur:
# 1) Plotting data asli
# 2) Low-pass filter
# 3) Offset waktu
# 4) Offset tinggi muka laut
# 5) Analisis harmonik (opsi 9 konstanta / konstanta lengkap)
# 6) Penghitungan datum pasut
# 7) Penghitungan tren kenaikan muka laut
# 8) Auto-pilih sensor terbaik dari semua sensor muka laut
# 9) QC flag per data (completeness, outlier, gap, flatline)
# =============================================================

suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(lubridate)
  library(ggplot2)
  library(tidyr)
  library(stringr)
})

# -----------------------------
# 0) PARAMETER UTAMA
# -----------------------------
CFG <- list(
  file_path = "data_pasut.xlsx",            # path data input
  sheet = 1,                                 # sheet excel
  timestamp_col = "Timestamp",              # kolom waktu
  input_tz = "UTC",                         # timezone input

  # Sensor
  auto_select_sensor = TRUE,                 # TRUE = otomatis pilih sensor terbaik
  sea_level_col = "PRS1 (m)",               # dipakai bila auto_select_sensor = FALSE
  sensor_pattern = "(PRS|RAD).*(m)",        # pola kandidat sensor muka laut

  # Offset
  time_offset_hours = 0,
  level_offset_m = 0.0,

  # QC
  outlier_method = "hampel",                # hampel / zscore
  hampel_window_n = 31,                      # jumlah titik (ganjil)
  hampel_n_sigma = 3,
  zscore_threshold = 3,
  flatline_tol = 1e-4,                       # toleransi flatline
  max_gap_multiplier = 3,                    # gap bila delta_t > multiplier * dt_median

  # Low-pass
  lowpass_method = "moving_average",        # moving_average / butterworth
  ma_window_hours = 25,
  butter_order = 4,
  butter_cutoff_hours = 30,

  # Harmonik: "9" atau "full"
  harmonic_constituents_mode = "9",

  out_dir = "output_pasut"
)

# -----------------------------
# 1) UTILITAS
# -----------------------------
ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

read_tide_data <- function(file_path, sheet = 1) {
  ext <- tolower(tools::file_ext(file_path))
  if (ext %in% c("xlsx", "xls")) {
    df <- readxl::read_excel(file_path, sheet = sheet)
  } else if (ext %in% c("csv", "txt")) {
    df <- readr::read_csv(file_path, show_col_types = FALSE)
  } else {
    stop("Format file belum didukung. Gunakan xlsx/xls/csv/txt.")
  }
  as_tibble(df)
}

infer_dt_minutes <- function(time_vec) {
  dt <- median(as.numeric(diff(time_vec), units = "mins"), na.rm = TRUE)
  if (is.na(dt) || dt <= 0) stop("Gagal mengestimasi interval waktu data.")
  dt
}

hampel_outlier <- function(x, k = 15, n_sigma = 3) {
  n <- length(x)
  flags <- rep(FALSE, n)
  if (n < (2 * k + 1)) return(flags)

  for (i in seq_len(n)) {
    left <- max(1, i - k)
    right <- min(n, i + k)
    w <- x[left:right]
    med <- median(w, na.rm = TRUE)
    madv <- mad(w, center = med, constant = 1, na.rm = TRUE)
    if (!is.na(x[i]) && !is.na(madv) && madv > 0) {
      flags[i] <- abs(x[i] - med) > (n_sigma * 1.4826 * madv)
    }
  }
  flags
}

zscore_outlier <- function(x, threshold = 3) {
  mu <- mean(x, na.rm = TRUE)
  sdv <- sd(x, na.rm = TRUE)
  if (is.na(sdv) || sdv == 0) return(rep(FALSE, length(x)))
  abs((x - mu) / sdv) > threshold
}

make_qc_flags <- function(time_vec, x, cfg) {
  dt_median <- infer_dt_minutes(time_vec)
  dt_vec <- c(NA_real_, as.numeric(diff(time_vec), units = "mins"))

  flag_missing <- is.na(x)

  if (cfg$outlier_method == "hampel") {
    k <- max(3, floor(cfg$hampel_window_n / 2))
    flag_outlier <- hampel_outlier(x, k = k, n_sigma = cfg$hampel_n_sigma)
  } else {
    flag_outlier <- zscore_outlier(x, threshold = cfg$zscore_threshold)
  }
  flag_outlier[is.na(flag_outlier)] <- FALSE

  dx <- c(NA_real_, abs(diff(x)))
  flag_flatline <- !is.na(dx) & dx <= cfg$flatline_tol

  flag_gap <- !is.na(dt_vec) & dt_vec > (cfg$max_gap_multiplier * dt_median)

  qc_code <- case_when(
    flag_missing ~ "MISSING",
    flag_outlier ~ "OUTLIER",
    flag_gap ~ "GAP",
    flag_flatline ~ "FLATLINE",
    TRUE ~ "GOOD"
  )

  tibble(
    flag_missing = flag_missing,
    flag_outlier = flag_outlier,
    flag_gap = flag_gap,
    flag_flatline = flag_flatline,
    qc_code = qc_code
  )
}

summarize_sensor_quality <- function(time_vec, x, cfg, sensor_name) {
  qc <- make_qc_flags(time_vec, x, cfg)
  n_total <- length(x)
  n_non_na <- sum(!qc$flag_missing)
  completeness <- ifelse(n_total > 0, n_non_na / n_total, 0)
  n_outlier <- sum(qc$flag_outlier, na.rm = TRUE)
  n_gap <- sum(qc$flag_gap, na.rm = TRUE)
  n_flatline <- sum(qc$flag_flatline, na.rm = TRUE)

  # Skor lebih kecil = lebih baik
  # Prioritas user: outlier dan completeness
  score <- (n_outlier / max(n_non_na, 1)) + (1 - completeness)

  tibble(
    sensor = sensor_name,
    n_total = n_total,
    n_non_na = n_non_na,
    completeness = completeness,
    n_outlier = n_outlier,
    n_gap = n_gap,
    n_flatline = n_flatline,
    score = score
  )
}

select_best_sensor <- function(df, time_col, cfg) {
  sensor_candidates <- names(df)[str_detect(names(df), cfg$sensor_pattern)]
  sensor_candidates <- sensor_candidates[sensor_candidates != time_col]

  if (length(sensor_candidates) == 0) {
    stop("Tidak ada sensor kandidat yang cocok dengan pola sensor_pattern.")
  }

  time_vec <- df[[time_col]]
  qlist <- lapply(sensor_candidates, function(s) {
    x <- suppressWarnings(as.numeric(df[[s]]))
    summarize_sensor_quality(time_vec, x, cfg, s)
  })

  quality <- bind_rows(qlist) |>
    arrange(score, desc(completeness), n_outlier)

  best <- quality$sensor[1]
  list(best_sensor = best, quality_table = quality)
}

lowpass_moving_average <- function(x, window_n) {
  stats::filter(x, rep(1 / window_n, window_n), sides = 2) |> as.numeric()
}

lowpass_butterworth <- function(x, dt_minutes, order = 4, cutoff_hours = 30) {
  if (!requireNamespace("signal", quietly = TRUE)) {
    stop("Paket 'signal' belum terpasang. Install: install.packages('signal')")
  }
  fs <- 60 / dt_minutes
  fc <- 1 / cutoff_hours
  wn <- fc / (fs / 2)
  if (wn >= 1) stop("Cutoff terlalu tinggi untuk sampling saat ini.")

  bf <- signal::butter(order, wn, type = "low")
  signal::filtfilt(bf, x)
}

find_tidal_extremes <- function(df, col = "sea_level_corr") {
  y <- df[[col]]
  idx <- which(!is.na(y))
  if (length(idx) < 5) return(list(high = integer(0), low = integer(0)))

  dy <- diff(y)
  sign_change <- diff(sign(dy))
  high_idx <- which(sign_change < 0) + 1
  low_idx <- which(sign_change > 0) + 1

  list(high = high_idx, low = low_idx)
}

compute_datums <- function(df, col = "sea_level_corr") {
  ex <- find_tidal_extremes(df, col)
  high_vals <- df[[col]][ex$high]
  low_vals <- df[[col]][ex$low]
  all_vals <- df[[col]]

  mhw <- mean(high_vals, na.rm = TRUE)
  mlw <- mean(low_vals, na.rm = TRUE)
  msl <- mean(all_vals, na.rm = TRUE)
  mtl <- (mhw + mlw) / 2
  mtr <- mhw - mlw
  hhat <- max(high_vals, na.rm = TRUE)
  llat <- min(low_vals, na.rm = TRUE)

  tibble(
    datum = c("MSL", "MHW", "MLW", "MTL", "MTR", "HHAT", "LLAT"),
    value_m = c(msl, mhw, mlw, mtl, mtr, hhat, llat)
  )
}

build_constituent_set <- function(mode = "9") {
  mode <- tolower(mode)
  if (mode == "9") {
    # 9 konstanta umum
    return(c("M2", "S2", "N2", "K2", "K1", "O1", "P1", "Q1", "M4"))
  }
  # "full" -> serahkan ke default TideHarmonics (fit sebanyak mungkin)
  NULL
}

harmonic_analysis <- function(df, time_col = "time_corr", level_col = "sea_level_corr", mode = "9") {
  if (!requireNamespace("TideHarmonics", quietly = TRUE)) {
    warning("Paket 'TideHarmonics' belum terpasang. Lewati analisis harmonik.\nInstall: install.packages('TideHarmonics')")
    return(NULL)
  }

  d2 <- df |>
    select(time = all_of(time_col), h = all_of(level_col), qc_code) |>
    filter(qc_code != "MISSING") |>
    tidyr::drop_na(time, h)

  if (nrow(d2) < 200) {
    warning("Data terlalu pendek untuk analisis harmonik yang stabil.")
  }

  hcn <- build_constituent_set(mode)
  if (is.null(hcn)) {
    fit <- TideHarmonics::ftide(x = d2$h, dto = d2$time)
  } else {
    fit <- TideHarmonics::ftide(x = d2$h, dto = d2$time, hcn = hcn)
  }

  coef <- as.data.frame(fit$coef)
  coef$constituent <- rownames(coef)
  rownames(coef) <- NULL
  list(fit = fit, coef = as_tibble(coef), n_used = nrow(d2), mode = mode)
}

plot_timeseries <- function(df, out_dir, sensor_name) {
  p1 <- ggplot(df, aes(time_corr, sea_level_raw)) +
    geom_line(alpha = 0.7) +
    geom_point(data = df |> filter(qc_code != "GOOD"), aes(time_corr, sea_level_raw, color = qc_code), size = 0.8) +
    labs(
      title = sprintf("Tinggi Muka Laut Raw - Sensor: %s", sensor_name),
      x = "Waktu", y = "Elevasi (m)", color = "QC"
    ) +
    theme_minimal(base_size = 12)

  p2 <- ggplot(df, aes(time_corr)) +
    geom_line(aes(y = sea_level_corr, color = "Offset"), alpha = 0.8) +
    geom_line(aes(y = sea_level_lp, color = "Low-pass"), linewidth = 0.7) +
    labs(title = "Raw (sudah offset) vs Low-pass", x = "Waktu", y = "Elevasi (m)", color = "Seri") +
    theme_minimal(base_size = 12)

  ggsave(file.path(out_dir, "01_raw_qc_timeseries.png"), p1, width = 11, height = 4, dpi = 150)
  ggsave(file.path(out_dir, "02_lowpass_timeseries.png"), p2, width = 11, height = 4, dpi = 150)
}

sea_level_trend <- function(df, time_col = "time_corr", level_col = "sea_level_lp") {
  d <- df |>
    filter(qc_code != "MISSING") |>
    select(time = all_of(time_col), y = all_of(level_col)) |>
    drop_na()

  monthly <- d |>
    mutate(month = floor_date(time, unit = "month")) |>
    group_by(month) |>
    summarise(msl_month = mean(y, na.rm = TRUE), .groups = "drop")

  t_year <- as.numeric(difftime(monthly$month, min(monthly$month), units = "days")) / 365.25
  fit <- lm(msl_month ~ t_year, data = mutate(monthly, t_year = t_year))

  slope_m_per_year <- coef(fit)[["t_year"]]
  slope_mm_per_year <- slope_m_per_year * 1000

  list(monthly = monthly, fit = fit, slope_mm_per_year = slope_mm_per_year)
}

plot_trend <- function(tr, out_dir) {
  p <- ggplot(tr$monthly, aes(month, msl_month)) +
    geom_point(size = 1.2, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "red") +
    labs(
      title = sprintf("Tren Kenaikan Muka Laut (%.2f mm/tahun)", tr$slope_mm_per_year),
      x = "Waktu", y = "MSL bulanan (m)"
    ) +
    theme_minimal(base_size = 12)

  ggsave(file.path(out_dir, "03_sea_level_trend.png"), p, width = 11, height = 4, dpi = 150)
}

# -----------------------------
# 2) ALUR UTAMA
# -----------------------------
main <- function(cfg) {
  ensure_dir(cfg$out_dir)

  message("[1/9] Membaca data...")
  df <- read_tide_data(cfg$file_path, cfg$sheet)
  stopifnot(cfg$timestamp_col %in% names(df))

  message("[2/9] Parsing waktu...")
  df <- df |>
    mutate(
      time_raw = parse_date_time(
        .data[[cfg$timestamp_col]],
        orders = c("dmy HM", "dmy HMS", "ymd HMS", "ymd HM", "mdy HM", "mdy HMS"),
        tz = cfg$input_tz
      )
    ) |>
    arrange(time_raw)

  message("[3/9] Seleksi sensor terbaik...")
  if (isTRUE(cfg$auto_select_sensor)) {
    sel <- select_best_sensor(df, "time_raw", cfg)
    selected_sensor <- sel$best_sensor
    sensor_quality <- sel$quality_table
  } else {
    if (!(cfg$sea_level_col %in% names(df))) {
      stop("sea_level_col tidak ditemukan pada data input.")
    }
    selected_sensor <- cfg$sea_level_col
    sensor_quality <- summarize_sensor_quality(
      time_vec = df$time_raw,
      x = suppressWarnings(as.numeric(df[[selected_sensor]])),
      cfg = cfg,
      sensor_name = selected_sensor
    )
  }

  message(sprintf("Sensor terpilih: %s", selected_sensor))
  readr::write_csv(sensor_quality, file.path(cfg$out_dir, "00_sensor_quality_ranking.csv"))

  message("[4/9] Bangun seri muka laut + offset...")
  df <- df |>
    mutate(
      sea_level_raw = as.numeric(.data[[selected_sensor]]),
      time_corr = time_raw + hours(cfg$time_offset_hours),
      sea_level_corr = sea_level_raw + cfg$level_offset_m,
      selected_sensor = selected_sensor
    ) |>
    arrange(time_corr)

  message("[5/9] QC flag per data...")
  qc_tbl <- make_qc_flags(df$time_corr, df$sea_level_corr, cfg)
  df <- bind_cols(df, qc_tbl)
  qc_summary <- df |>
    count(qc_code, name = "n") |>
    mutate(pct = n / sum(n) * 100)
  readr::write_csv(qc_summary, file.path(cfg$out_dir, "01_qc_summary.csv"))

  dt_minutes <- infer_dt_minutes(df$time_corr)
  message(sprintf("Interval sampling terdeteksi: %.2f menit", dt_minutes))

  message("[6/9] Low-pass filtering...")
  x_for_filter <- ifelse(df$qc_code %in% c("MISSING", "OUTLIER"), NA_real_, df$sea_level_corr)

  if (cfg$lowpass_method == "moving_average") {
    window_n <- max(3, round(cfg$ma_window_hours * 60 / dt_minutes))
    if (window_n %% 2 == 0) window_n <- window_n + 1
    df$sea_level_lp <- lowpass_moving_average(x_for_filter, window_n)
  } else if (cfg$lowpass_method == "butterworth") {
    # interpolasi sederhana agar filtfilt stabil terhadap NA
    idx_valid <- which(!is.na(x_for_filter))
    if (length(idx_valid) < 3) {
      stop("Data valid terlalu sedikit untuk Butterworth filtering.")
    }
    x_interp <- approx(
      x = idx_valid,
      y = x_for_filter[idx_valid],
      xout = seq_along(x_for_filter),
      method = "linear",
      rule = 2
    )$y
    df$sea_level_lp <- lowpass_butterworth(
      x = x_interp,
      dt_minutes = dt_minutes,
      order = cfg$butter_order,
      cutoff_hours = cfg$butter_cutoff_hours
    )
  } else {
    stop("Metode lowpass tidak dikenali. Gunakan 'moving_average' atau 'butterworth'.")
  }

  message("[7/9] Plotting...")
  plot_timeseries(df, cfg$out_dir, selected_sensor)

  message("[8/9] Analisis harmonik + datum...")
  har <- harmonic_analysis(df, "time_corr", "sea_level_corr", cfg$harmonic_constituents_mode)
  if (!is.null(har)) {
    readr::write_csv(har$coef, file.path(cfg$out_dir, "04_harmonic_constituents.csv"))
    readr::write_csv(
      tibble(metric = c("mode", "n_used"), value = c(har$mode, har$n_used)),
      file.path(cfg$out_dir, "04_harmonic_metadata.csv")
    )
  }

  datums <- compute_datums(df, "sea_level_corr")
  readr::write_csv(datums, file.path(cfg$out_dir, "05_tidal_datums.csv"))

  message("[9/9] Tren kenaikan muka laut...")
  tr <- sea_level_trend(df, "time_corr", "sea_level_lp")
  plot_trend(tr, cfg$out_dir)

  trend_summary <- tibble(
    metric = c("selected_sensor", "slope_mm_per_year", "n_months"),
    value = c(selected_sensor, as.character(tr$slope_mm_per_year), as.character(nrow(tr$monthly)))
  )
  readr::write_csv(trend_summary, file.path(cfg$out_dir, "06_sea_level_trend_summary.csv"))

  readr::write_csv(df, file.path(cfg$out_dir, "00_processed_timeseries.csv"))

  message("Selesai. Output tersimpan di folder: ", cfg$out_dir)
}

if (sys.nframe() == 0) {
  # Contoh:
  # Rscript scripts_analisis_pasut.R data.xlsx
  # Rscript scripts_analisis_pasut.R data.xlsx "full"
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) >= 1) CFG$file_path <- args[[1]]
  if (length(args) >= 2) CFG$harmonic_constituents_mode <- args[[2]]

  main(CFG)
}
