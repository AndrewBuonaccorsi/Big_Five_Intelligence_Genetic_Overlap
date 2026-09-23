#!/usr/bin/env Rscript


`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

DEFAULTS <- list(
  base            = ".",
  cache_dir       = "",
  out_dir         = "",
  full_subdir     = "geneset_output",
  traits          = "iq,extra,agree,consc,neurot,open",
  ranking_file    = "",
  fdr             = "0.05",
  fdr_scope       = "trait",          # trait | trait_collection
  n_label         = "10",
  label_maxchar   = "45",
  label_cex       = "0.62",
  collections     = "all",
  skip_figures    = "FALSE",
  no_ngenes       = "FALSE"
)

parse_args <- function(argv) {
  opt <- DEFAULTS; i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a, call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", a), fixed = TRUE)
    if (!key %in% names(DEFAULTS)) stop("Unknown option: --", sub("^--", "", a), call. = FALSE)
    is_flag <- key %in% c("skip_figures", "no_ngenes")
    nxt <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
    if (is_flag && (is.na(nxt) || startsWith(nxt, "--"))) { opt[[key]] <- "TRUE"; i <- i + 1L }
    else {
      if (is.na(nxt)) stop("Option --", sub("^--", "", a), " needs a value", call. = FALSE)
      opt[[key]] <- nxt; i <- i + 2L
    }
  }
  opt
}

opt      <- parse_args(commandArgs(trailingOnly = TRUE))
BASE     <- normalizePath(opt$base, mustWork = TRUE)
CACHEDIR <- opt$cache_dir %||% file.path(BASE, "results", "cache")
OUTDIR   <- opt$out_dir   %||% file.path(BASE, "results", "enrich_diff")
FULLROOT <- file.path(BASE, opt$full_subdir)
TRAITS   <- trimws(strsplit(opt$traits, ",", fixed = TRUE)[[1]])
NT       <- length(TRAITS)
FDR      <- as.numeric(opt$fdr)
SCOPE    <- tolower(opt$fdr_scope)
NLAB     <- as.integer(opt$n_label)
LABMAX   <- as.integer(opt$label_maxchar)
LABCEX   <- as.numeric(opt$label_cex)
SKIPFIG  <- toupper(opt$skip_figures) %in% c("TRUE", "T", "1", "YES")
NONGENES <- toupper(opt$no_ngenes)    %in% c("TRUE", "T", "1", "YES")
RANKFILE <- opt$ranking_file %||% file.path(BASE, "results", "gene_set_ranking.txt")

if (!SCOPE %in% c("trait", "trait_collection"))
  stop("--fdr-scope must be 'trait' or 'trait_collection'", call. = FALSE)
if (NT != 6L)
  stop("This analysis is defined for 6 traits; got ", NT, ".", call. = FALSE)

FIGDIR <- file.path(OUTDIR, "figures")
ALLDIR <- file.path(OUTDIR, "all_results")
SIGDIR <- file.path(OUTDIR, "significant_results")
for (d in c(OUTDIR, FIGDIR, ALLDIR, SIGDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

TRAIT_LABELS <- c(iq = "Intelligence", extra = "Extraversion", agree = "Agreeableness",
                  consc = "Conscientiousness", neurot = "Neuroticism", open = "Openness")
PANEL_ORDER <- intersect(names(TRAIT_LABELS), TRAITS)
if (length(PANEL_ORDER) != NT) PANEL_ORDER <- TRAITS
lab_of <- function(tr) if (!is.na(TRAIT_LABELS[tr])) unname(TRAIT_LABELS[tr]) else tr

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = ""); flush.console()
}
DIAG <- new.env(parent = emptyenv()); DIAG$lines <- character(0)
diag_add <- function(...) DIAG$lines <- c(DIAG$lines, paste0(...))

## ------------------------------- read the cache ------------------------------

metas <- list.files(CACHEDIR, pattern = "\\.meta\\.rds$", full.names = TRUE)
if (!length(metas))
  stop("No cached collections in ", CACHEDIR, ".\n",
       "  Build the cache first:  ./run_ranker.sh --stage parse --resume", call. = FALSE)

meta  <- lapply(metas, readRDS)
colls <- vapply(meta, `[[`, character(1), "collection")
ord   <- order(colls); meta <- meta[ord]; colls <- colls[ord]

if (!identical(opt$collections, "all")) {
  want  <- trimws(strsplit(opt$collections, ",", fixed = TRUE)[[1]])
  keep  <- colls %in% want
  if (!any(keep)) stop("None of the requested collections are cached.", call. = FALSE)
  meta <- meta[keep]; colls <- colls[keep]
}

reps_ok <- Reduce(`&`, lapply(meta, `[[`, "reps_ok"))
NJK     <- length(reps_ok)
B       <- sum(reps_ok)
if (B < 2L) stop("Only ", B, " replicate(s) usable in every collection.", call. = FALSE)
slices  <- c(1L, 1L + which(reps_ok))       # slice 1 = full data
S       <- length(slices)
N       <- sum(vapply(meta, `[[`, numeric(1), "n_sets"))

say("cache                : ", CACHEDIR)
say("collections pooled   : ", length(colls))
say("gene sets pooled     : ", N)
say("replicates usable    : ", B, " of ", NJK)
say("traits (panel order) : ", paste(PANEL_ORDER, collapse = ", "))
say("FDR                  : ", FDR, "  (BH within ",
    if (SCOPE == "trait") "focal trait, all collections pooled" else "focal trait x collection", ")")

beta_cache <- function(coll) file.path(CACHEDIR, paste0(coll, ".beta.rds"))
read_beta  <- function(coll) {
  o <- readRDS(beta_cache(coll))
  if (!identical(sort(o$traits), sort(TRAITS)))
    stop("Cached traits for ", coll, " differ from --traits.", call. = FALSE)
  o
}
TR_IDX <- NULL


Ssum <- matrix(0, nrow = NT, ncol = S)
for (cc in colls) {
  o <- read_beta(cc)
  if (is.null(TR_IDX)) TR_IDX <- match(TRAITS, o$traits)
  Ssum <- Ssum + colSums(o$beta[, TR_IDX, slices, drop = FALSE], dims = 1L)
  rm(o)
}
mu <- Ssum / N

Qsum <- matrix(0, nrow = NT, ncol = S)
for (cc in colls) {
  o <- read_beta(cc)
  d <- sweep(o$beta[, TR_IDX, slices, drop = FALSE], c(2L, 3L), mu, "-")
  Qsum <- Qsum + colSums(d^2, dims = 1L)
  rm(o, d)
}
sdev <- sqrt(Qsum / (N - 1))
if (any(!is.finite(sdev)) || any(sdev == 0))
  stop("Degenerate pooled SD for some trait/replicate.", call. = FALSE)

diag_add("=== POOLED NORMALIZATION (full data) ===")
for (i in seq_len(NT))
  diag_add(sprintf("  %-8s mean BETA_STD = %12.6g   sd = %12.6g",
                   TRAITS[i], mu[i, 1L], sdev[i, 1L]))

parts <- vector("list", length(colls))
meanz <- vector("list", length(colls))

for (j in seq_along(colls)) {
  cc <- colls[[j]]
  o  <- read_beta(cc)
  n  <- length(o$sets)
  
  z <- sweep(sweep(o$beta[, TR_IDX, slices, drop = FALSE], c(2L, 3L), mu, "-"),
             c(2L, 3L), sdev, "/")          
  
  zsum <- z[, 1L, ]
  for (i in 2:NT) zsum <- zsum + z[, i, ]  
  zbar <- zsum / NT
  
  meanz[[j]] <- data.frame(COLLECTION = cc, GENE_SET = o$sets,
                           MEAN_Z = zbar[, 1L], stringsAsFactors = FALSE)
  
  rows <- vector("list", NT)
  for (i in seq_len(NT)) {
    zi   <- z[, i, ]                       
    comp <- (zsum - zi) / (NT - 1L)               
    D    <- zi - comp                          
    
    Dfull <- D[, 1L]
    Dj    <- D[, -1L, drop = FALSE]
    Dbar  <- rowMeans(Dj)
    seD   <- sqrt(((B - 1) / B) * rowSums((Dj - Dbar)^2))
    
    rows[[i]] <- data.frame(
      COLLECTION     = cc,
      GENE_SET       = o$sets,
      TRAIT          = TRAITS[i],
      TRAIT_LABEL    = lab_of(TRAITS[i]),
      Z_FOCAL        = zi[, 1L],
      MEAN_Z_OTHERS  = comp[, 1L],
      MEAN_Z_ALL6    = zbar[, 1L],
      D              = Dfull,
      SE_D           = seD,
      stringsAsFactors = FALSE
    )
  }
  parts[[j]] <- do.call(rbind, rows)
  rm(o, z, zsum, zbar, rows); invisible(gc(FALSE))
  say(sprintf("  %-42s  sets=%6d", cc, n))
}

res <- do.call(rbind, parts); rm(parts)
meanz <- do.call(rbind, meanz)

res$STAT <- res$D / res$SE_D
res$P    <- 2 * stats::pnorm(-abs(res$STAT))

bad <- !is.finite(res$STAT)
if (any(bad)) {
  diag_add("", "  WARN ", sum(bad), " row(s) with non-finite STAT (zero or NA SE); P set to NA")
  res$P[bad] <- NA_real_
}


fam <- if (SCOPE == "trait") res$TRAIT else paste(res$TRAIT, res$COLLECTION, sep = "\u001f")
res$P_BH <- NA_real_
for (f in unique(fam)) {
  k <- which(fam == f)
  res$P_BH[k] <- stats::p.adjust(res$P[k], method = "BH")
}
res$SIG <- !is.na(res$P_BH) & res$P_BH < FDR


if (file.exists(RANKFILE)) {
  rk <- utils::read.delim(RANKFILE, stringsAsFactors = FALSE, quote = "", comment.char = "")
  if (all(c("COLLECTION", "GENE_SET", "MEAN_Z") %in% names(rk))) {
    k1 <- paste(meanz$COLLECTION, meanz$GENE_SET, sep = "\u001f")
    k2 <- paste(rk$COLLECTION,    rk$GENE_SET,    sep = "\u001f")
    m  <- match(k1, k2)
    ok <- !is.na(m)
    if (any(ok)) {
      dmax <- max(abs(meanz$MEAN_Z[ok] - as.numeric(rk$MEAN_Z[m[ok]])))
      diag_add("", "=== CONSISTENCY WITH gene_set_ranking.txt ===")
      diag_add("  gene sets matched      : ", sum(ok), " of ", nrow(meanz))
      diag_add("  max |MEAN_Z difference| : ", sprintf("%.3e", dmax),
               if (dmax < 1e-5) "   (rounding only -- identical z scores)" else
                 "   *** MISMATCH: the cache and the ranking disagree ***")
      say("consistency vs ranking: max |dMEAN_Z| = ", sprintf("%.3e", dmax),
          " over ", sum(ok), " matched gene sets")
      if (dmax >= 1e-5)
        say("WARNING: MEAN_Z does not reproduce gene_set_ranking.txt -- ",
            "the cache may have changed since that file was written.")
    }
  }
} else {
  diag_add("", "  NOTE ", RANKFILE, " not found; skipped the consistency check.")
}


if (!NONGENES) {
  ng <- tryCatch({
    read_ngenes <- function(path) {
      ln <- readLines(path, warn = FALSE)
      ln <- trimws(ln[!grepl("^\\s*#", ln)]); ln <- ln[nzchar(ln)]
      h  <- grep("^VARIABLE\\b", ln)[1]
      if (is.na(h) || h >= length(ln)) return(NULL)
      hdr <- strsplit(ln[h], "[ \t]+")[[1]]; k <- length(hdr)
      nc <- match("FULL_NAME", hdr); if (is.na(nc)) nc <- match("VARIABLE", hdr)
      gc_ <- match("NGENES", hdr); tc <- match("TYPE", hdr)
      if (is.na(gc_)) return(NULL)
      sp <- strsplit(ln[(h + 1L):length(ln)], "[ \t]+")
      sp <- sp[lengths(sp) >= k]
      if (!length(sp)) return(NULL)
      m <- if (all(lengths(sp) == k)) matrix(unlist(sp, use.names = FALSE), nrow = k) else
        vapply(sp, function(x) if (length(x) > k)
          c(x[seq_len(k - 1L)], paste(x[k:length(x)], collapse = " ")) else x, character(k))
      keep <- if (!is.na(tc)) m[tc, ] == "SET" else rep(TRUE, ncol(m))
      data.frame(GENE_SET = m[nc, keep],
                 NGENES = suppressWarnings(as.numeric(m[gc_, keep])),
                 stringsAsFactors = FALSE)
    }
    out <- list()
    for (cc in colls) for (tr in TRAITS) {
      p <- file.path(FULLROOT, cc, paste0(tr, ".gsa.out"))
      if (!file.exists(p)) next
      d <- read_ngenes(p); if (is.null(d)) next
      d <- d[!duplicated(d$GENE_SET), , drop = FALSE]
      d$COLLECTION <- cc; d$TRAIT <- tr
      out[[length(out) + 1L]] <- d
    }
    if (!length(out)) NULL else do.call(rbind, out)
  }, error = function(e) NULL)
  
  if (!is.null(ng)) {
    k1 <- paste(res$COLLECTION, res$GENE_SET, res$TRAIT, sep = "\u001f")
    k2 <- paste(ng$COLLECTION,  ng$GENE_SET,  ng$TRAIT,  sep = "\u001f")
    res$NGENES <- ng$NGENES[match(k1, k2)]
    diag_add("", "  NGENES attached for ",
             sprintf("%.1f%%", 100 * mean(!is.na(res$NGENES))), " of rows")
  } else {
    res$NGENES <- NA_real_
    diag_add("", "  NGENES could not be read; column left as NA")
  }
} else res$NGENES <- NA_real_


COLS <- c("COLLECTION", "GENE_SET", "TRAIT", "TRAIT_LABEL", "NGENES",
          "Z_FOCAL", "MEAN_Z_OTHERS", "MEAN_Z_ALL6", "D", "SE_D", "STAT",
          "P", "P_BH", "SIG")
res <- res[, COLS, drop = FALSE]

utils::write.table(res[order(res$TRAIT, -abs(res$STAT)), , drop = FALSE],
                   file.path(OUTDIR, "all_identity_deviations.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)

sig_all <- res[res$SIG, , drop = FALSE]
sig_all <- sig_all[order(sig_all$COLLECTION, sig_all$TRAIT, -abs(sig_all$D)), , drop = FALSE]
utils::write.csv(sig_all,
                 file.path(SIGDIR, "ALL_COLLECTIONS_SIGNIFICANT_identity_deviations.csv"),
                 row.names = FALSE, na = "NA")

for (cc in colls) {
  d <- res[res$COLLECTION == cc, , drop = FALSE]
  d <- d[order(d$TRAIT, -abs(d$D)), , drop = FALSE]
  utils::write.csv(d, file.path(ALLDIR, paste0(cc, "_ALL_identity_deviations.csv")),
                   row.names = FALSE, na = "NA")
  s <- d[d$SIG, , drop = FALSE]
  utils::write.csv(s, file.path(SIGDIR, paste0(cc, "_SIGNIFICANT_identity_deviations.csv")),
                   row.names = FALSE, na = "NA")
}

place_labels <- function(x, y, labs, cex = LABCEX) {
  if (!length(x)) return(invisible(NULL))
  usr <- par("usr"); pin <- par("pin")
  sx <- (usr[2] - usr[1]) / pin[1]; sy <- (usr[4] - usr[3]) / pin[2]
  w <- strwidth(labs, cex = cex); h <- strheight(labs, cex = cex)
  dirs <- seq(0, 2 * pi, length.out = 25)[-25]
  rads <- c(0.20, 0.32, 0.46, 0.62)
  placed <- list()
  for (k in order(-abs(y - x))) {
    best <- NULL; bestpen <- Inf
    for (rr in rads) for (a in dirs) {
      cx <- x[k] + cos(a) * rr * sx; cy <- y[k] + sin(a) * rr * sy
      r  <- c(cx - w[k] / 2, cx + w[k] / 2, cy - h[k] / 2, cy + h[k] / 2)
      pen <- 40 * ((max(0, usr[1] - r[1]) + max(0, r[2] - usr[2])) / sx +
                     (max(0, usr[3] - r[3]) + max(0, r[4] - usr[4])) / sy)
      for (p in placed) {
        ox <- max(0, min(r[2], p[2]) - max(r[1], p[1])) / sx
        oy <- max(0, min(r[4], p[4]) - max(r[3], p[3])) / sy
        pen <- pen + 60 * ox * oy
      }
      for (m in seq_along(x)) {
        if (x[m] >= r[1] && x[m] <= r[2] && y[m] >= r[3] && y[m] <= r[4]) pen <- pen + 3
      }
      pen <- pen + 0.6 * rr
      if (pen < bestpen) { bestpen <- pen; best <- r }
    }
    placed[[length(placed) + 1L]] <- best
    lx <- mean(best[1:2]); ly <- mean(best[3:4])
    segments(x[k], y[k], lx, ly, col = "grey45", lwd = 0.4)
    text(lx, ly, labs[k], cex = cex, col = "black", xpd = NA)
  }
  invisible(NULL)
}

shorten <- function(s, n) ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "..."), s)

six_panel <- function(d, path, main_title) {
  grDevices::pdf(path, width = 16, height = 10, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  par(mfrow = c(2, 3), mar = c(4.4, 4.4, 3.2, 1.4), oma = c(0, 0, 2.4, 0),
      mgp = c(2.5, 0.7, 0), las = 1, bty = "l")
  
  for (pi_ in seq_along(PANEL_ORDER)) {
    tr <- PANEL_ORDER[[pi_]]
    dd <- d[d$TRAIT == tr, , drop = FALSE]
    if (!nrow(dd)) { plot.new(); title(main = lab_of(tr)); next }
    
    xr <- range(dd$MEAN_Z_OTHERS, finite = TRUE)
    yr <- range(dd$Z_FOCAL, finite = TRUE)
    pad <- 0.10 * c(diff(xr), diff(yr))
    plot(dd$MEAN_Z_OTHERS, dd$Z_FOCAL, type = "n",
         xlim = xr + c(-1, 1) * pad[1], ylim = yr + c(-1, 1) * pad[2],
         xlab = "Mean enrichment of other traits",
         ylab = paste0(lab_of(tr), " enrichment"),
         main = lab_of(tr), cex.main = 1.25, font.main = 2, cex.lab = 1.0)
    grid(col = "grey92", lty = 1, lwd = 0.6)
    
    ns <- !dd$SIG
    points(dd$MEAN_Z_OTHERS[ns], dd$Z_FOCAL[ns], pch = 16, cex = 0.35,
           col = grDevices::adjustcolor("grey35", 0.30))
    abline(0, 1, lty = 2, lwd = 1.0, col = "grey20")
    fit <- stats::lm(Z_FOCAL ~ MEAN_Z_OTHERS, data = dd)
    abline(fit, col = "#2166AC", lwd = 1.3)
    if (any(dd$SIG)) {
      points(dd$MEAN_Z_OTHERS[dd$SIG], dd$Z_FOCAL[dd$SIG], pch = 21, cex = 0.72,
             bg = grDevices::adjustcolor("#B2182B", 0.90), col = "white", lwd = 0.5)
    }
    
    lb <- dd[dd$SIG, , drop = FALSE]
    if (nrow(lb)) {
      lb <- lb[order(-abs(lb$D)), , drop = FALSE]
      if (NLAB > 0L) lb <- lb[seq_len(min(NLAB, nrow(lb))), , drop = FALSE]
      place_labels(lb$MEAN_Z_OTHERS, lb$Z_FOCAL, shorten(lb$GENE_SET, LABMAX))
    }
    
    b <- unname(coef(fit)[2]); r <- stats::cor(dd$MEAN_Z_OTHERS, dd$Z_FOCAL)
    legend("bottomright", bty = "n", cex = 0.85,
           legend = c(sprintf("beta = %.2f   r = %.2f", b, r),
                      sprintf("%d of %d FDR < %.2f", sum(dd$SIG), nrow(dd), FDR)))
    mtext(LETTERS[pi_], side = 3, adj = -0.10, line = 1.1, font = 2, cex = 1.05)
  }
  mtext(main_title, side = 3, outer = TRUE, font = 2, cex = 1.15)
  invisible(NULL)
}

if (!SKIPFIG) {
  for (cc in colls) {
    six_panel(res[res$COLLECTION == cc, , drop = FALSE],
              file.path(FIGDIR, paste0(cc, "_six_panel_identity_plot.pdf")), cc)
    say("  figure: ", cc)
  }
  six_panel(res, file.path(FIGDIR, "ALL_COLLECTIONS_six_panel_identity_plot.pdf"),
            paste0("All collections (", format(N, big.mark = ","), " gene sets)"))
  say("  figure: ALL_COLLECTIONS")
} else say("figures skipped (--skip-figures)")


diag_add("", "=== SIGNIFICANT DIFFERENCES BY FOCAL TRAIT (BH < ", FDR, ") ===")
diag_add(sprintf("  %-20s %9s %9s %9s", "TRAIT", "N_TESTED", "N_SIG", "PCT"))
for (tr in PANEL_ORDER) {
  k <- res$TRAIT == tr
  diag_add(sprintf("  %-20s %9d %9d %8.2f%%", lab_of(tr), sum(k), sum(res$SIG[k]),
                   100 * sum(res$SIG[k]) / max(1, sum(k))))
}

diag_add("", "=== SIGNIFICANT DIFFERENCES BY COLLECTION ===")
diag_add(sprintf("  %-42s %9s %9s", "COLLECTION", "N_SIG", "N_TESTED"))
for (cc in colls) {
  k <- res$COLLECTION == cc
  diag_add(sprintf("  %-42s %9d %9d", cc, sum(res$SIG[k]), sum(k)))
}

diag_add("", "=== SUMMARY ===")
diag_add("  collections         : ", length(colls))
diag_add("  gene sets           : ", N)
diag_add("  tests               : ", nrow(res), "  (", N, " x 6 focal traits)")
diag_add("  replicates used (B) : ", B, " of ", NJK)
diag_add("  standardization     : pooled over all ", N, " gene sets in all collections")
diag_add("  BH family           : ",
         if (SCOPE == "trait") "focal trait (all collections pooled)" else "focal trait x collection")
diag_add("  significant at BH < ", FDR, " : ", sum(res$SIG))
diag_add("  NOTE the six panels sum to zero per gene set (5 df, not 6);")
diag_add("       the six tests for a gene set are not independent.")
diag_add("  finished            : ", format(Sys.time()))

writeLines(DIAG$lines, file.path(OUTDIR, "enrich_diff_diagnostics.txt"))

say("significant at BH < ", FDR, ": ", sum(res$SIG), " of ", nrow(res), " tests")
say("wrote ", OUTDIR)