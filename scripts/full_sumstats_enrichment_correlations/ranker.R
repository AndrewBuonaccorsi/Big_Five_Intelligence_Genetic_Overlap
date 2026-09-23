#!/usr/bin/env Rscript


`%||%` <- function(a, b) if (is.null(a) || is.na(a) || !nzchar(a)) b else a

DEFAULTS <- list(
  base                   = ".",
  jack_subdir            = "output",
  full_subdir            = "geneset_output",
  trait_dir_prefix       = "jack_",
  traits                 = "iq,open,consc,extra,agree,neurot",
  n_jk                   = "200",
  jk_pad                 = "3",
  collections            = "all",
  out                    = "results/gene_set_ranking.txt",
  cache_dir              = "",
  diagnostics            = "",
  stage                  = "all",
  rank_within_collection = "FALSE",
  resume                 = "FALSE",
  test                   = "FALSE"
)

parse_args <- function(argv) {
  opt <- DEFAULTS
  i <- 1L
  while (i <= length(argv)) {
    a <- argv[[i]]
    if (!startsWith(a, "--")) stop("Unexpected argument: ", a, call. = FALSE)
    key <- gsub("-", "_", sub("^--", "", a), fixed = TRUE)
    if (!key %in% names(DEFAULTS)) stop("Unknown option: --", sub("^--", "", a), call. = FALSE)
    is_flag <- key %in% c("rank_within_collection", "resume", "test")
    nxt <- if (i < length(argv)) argv[[i + 1L]] else NA_character_
    if (is_flag && (is.na(nxt) || startsWith(nxt, "--"))) {
      opt[[key]] <- "TRUE"; i <- i + 1L
    } else {
      if (is.na(nxt)) stop("Option --", sub("^--", "", a), " needs a value", call. = FALSE)
      opt[[key]] <- nxt; i <- i + 2L
    }
  }
  opt
}

opt      <- parse_args(commandArgs(trailingOnly = TRUE))
BASE     <- normalizePath(opt$base, mustWork = TRUE)
JACKROOT <- file.path(BASE, opt$jack_subdir)
FULLROOT <- file.path(BASE, opt$full_subdir)
TRAITS   <- trimws(strsplit(opt$traits, ",", fixed = TRUE)[[1]])
NT       <- length(TRAITS)
NJK      <- as.integer(opt$n_jk)
JKPAD    <- as.integer(opt$jk_pad)
OUTFILE  <- if (startsWith(opt$out, "/")) opt$out else file.path(BASE, opt$out)
OUTDIR   <- dirname(OUTFILE)
CACHEDIR <- opt$cache_dir %||% file.path(OUTDIR, "cache")
DIAGDIR  <- file.path(OUTDIR, "diagnostics")
STAGE    <- tolower(opt$stage)
RANK_WC  <- toupper(opt$rank_within_collection) %in% c("TRUE", "T", "1", "YES")
RESUME   <- toupper(opt$resume) %in% c("TRUE", "T", "1", "YES")
TESTMODE <- toupper(opt$test)   %in% c("TRUE", "T", "1", "YES")

if (!STAGE %in% c("all", "parse", "combine"))
  stop("--stage must be one of: all, parse, combine", call. = FALSE)

for (d in c(OUTDIR, CACHEDIR, DIAGDIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

say <- function(...) {
  cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")
  flush.console()
}

DIAG <- new.env(parent = emptyenv())
DIAG$lines <- character(0)
diag_add <- function(...) DIAG$lines <- c(DIAG$lines, paste0(...))


# Returns BETA_STD indexed by gene set name.
parse_gsa <- function(path) {
  ln <- tryCatch(readLines(path, warn = FALSE), error = function(e) NULL)
  if (is.null(ln) || !length(ln)) return(NULL)
  
  ln <- ln[!grepl("^\\s*#", ln)]
  ln <- trimws(ln)
  ln <- ln[nzchar(ln)]
  if (!length(ln)) return(NULL)
  
  hdr_i <- grep("^VARIABLE\\b", ln)
  if (!length(hdr_i)) return(NULL)
  hdr_i <- hdr_i[[1L]]
  if (hdr_i >= length(ln)) return(NULL)
  
  hdr <- strsplit(ln[hdr_i], "[ \t]+")[[1L]]
  k   <- length(hdr)
  
  beta_col <- match("BETA_STD", hdr)
  if (is.na(beta_col))
    stop("No BETA_STD column in ", path,
         " (header: ", paste(hdr, collapse = " "), ")", call. = FALSE)
  name_col <- match("FULL_NAME", hdr)
  name_src <- "FULL_NAME"
  if (is.na(name_col)) { name_col <- match("VARIABLE", hdr); name_src <- "VARIABLE" }
  type_col <- match("TYPE", hdr)
  
  dat  <- ln[(hdr_i + 1L):length(ln)]
  sp   <- strsplit(dat, "[ \t]+")
  lens <- lengths(sp)
  
  short <- lens < k
  if (any(short)) {
    diag_add("  WARN ", basename(path), ": dropped ", sum(short),
             " malformed row(s) with fewer than ", k, " fields")
    sp <- sp[!short]; lens <- lens[!short]
  }
  if (!length(sp)) return(NULL)
  
  if (all(lens == k)) {
    m <- matrix(unlist(sp, use.names = FALSE), nrow = k)
  } else {
    m <- vapply(sp, function(x) {
      if (length(x) > k) c(x[seq_len(k - 1L)], paste(x[k:length(x)], collapse = " ")) else x
    }, character(k))
  }
  
  keep <- rep(TRUE, ncol(m))
  if (!is.na(type_col)) keep <- m[type_col, ] == "SET"
  if (!any(keep)) return(NULL)
  
  nms  <- m[name_col, keep]
  beta <- suppressWarnings(as.numeric(m[beta_col, keep]))
  
  bad <- is.na(beta) | is.na(nms) | !nzchar(nms)
  if (any(bad)) {
    diag_add("  WARN ", basename(path), ": dropped ", sum(bad),
             " row(s) with unparseable BETA_STD or empty name")
    nms <- nms[!bad]; beta <- beta[!bad]
  }
  if (!length(nms)) return(NULL)
  
  dup <- duplicated(nms) | duplicated(nms, fromLast = TRUE)
  if (any(dup)) {
    diag_add("  WARN ", basename(path), ": dropped ", sum(dup),
             " row(s) with duplicated gene set name(s), e.g. ",
             paste(utils::head(unique(nms[dup]), 3L), collapse = ", "))
    nms <- nms[!dup]; beta <- beta[!dup]
  }
  
  v <- stats::setNames(beta, nms)
  attr(v, "name_col") <- name_src
  v
}

full_path <- function(coll, trait) file.path(FULLROOT, coll, paste0(trait, ".gsa.out"))
jk_path   <- function(coll, trait, b) {
  file.path(JACKROOT, coll, paste0(opt$trait_dir_prefix, trait),
            sprintf("%s_jk%0*d.gsa.out", trait, JKPAD, b))
}
beta_cache <- function(coll) file.path(CACHEDIR, paste0(coll, ".beta.rds"))
meta_cache <- function(coll) file.path(CACHEDIR, paste0(coll, ".meta.rds"))

list_subdirs <- function(p) {
  if (!dir.exists(p)) character(0) else
    basename(list.dirs(p, recursive = FALSE, full.names = TRUE))
}

if (identical(opt$collections, "all")) {
  jd <- list_subdirs(JACKROOT); fd <- list_subdirs(FULLROOT)
  COLLECTIONS <- sort(intersect(jd, fd))
  miss <- setdiff(union(jd, fd), COLLECTIONS)
  if (length(miss)) {
    say("NOTE: ", length(miss), " collection(s) present in only one of ",
        opt$jack_subdir, "/ and ", opt$full_subdir, "/ -- skipped: ",
        paste(miss, collapse = ", "))
    diag_add("SKIPPED (present in only one root): ", paste(miss, collapse = ", "))
  }
} else {
  COLLECTIONS <- trimws(strsplit(opt$collections, ",", fixed = TRUE)[[1]])
}
if (!length(COLLECTIONS)) stop("No collections found under ", JACKROOT, " and ", FULLROOT)

DIAGFILE <- if (nzchar(opt$diagnostics)) {
  opt$diagnostics
} else if (STAGE == "parse" && length(COLLECTIONS) == 1L) {
  file.path(DIAGDIR, paste0(COLLECTIONS[[1L]], ".txt"))
} else {
  file.path(OUTDIR, "gene_set_ranking_diagnostics.txt")
}


if (TESTMODE) {
  coll <- COLLECTIONS[[1L]]
  say("TEST MODE -- collection: ", coll)
  for (tr in TRAITS) {
    for (p in c(full_path(coll, tr), jk_path(coll, tr, 1L))) {
      if (!file.exists(p)) { cat("  MISSING  ", p, "\n"); next }
      v <- parse_gsa(p)
      if (is.null(v)) { cat("  UNPARSED ", p, "\n"); next }
      cat(sprintf("  %-58s n=%5d  names=%-9s  head: %s\n",
                  sub(paste0("^", BASE, "/"), "", p), length(v),
                  attr(v, "name_col"),
                  paste(sprintf("%s=%.5g", names(v)[1:min(2, length(v))],
                                unname(v)[1:min(2, length(v))]), collapse = "  ")))
    }
  }
  say("Collections detected (", length(COLLECTIONS), "): ",
      paste(COLLECTIONS, collapse = ", "))
  if (length(DIAG$lines)) cat(paste(DIAG$lines, collapse = "\n"), "\n", sep = "")
  say("Test complete. No results written.")
  quit(save = "no", status = 0L)
}



parse_collection <- function(coll) {
  t0 <- Sys.time()
  diag_add("", "=== ", coll, " ===")
  
  fpaths <- vapply(TRAITS, function(tr) full_path(coll, tr), character(1))
  if (!all(file.exists(fpaths))) {
    say("SKIP ", coll, " -- missing full-data file(s)")
    diag_add("  SKIPPED: missing full-data file(s): ",
             paste(basename(fpaths[!file.exists(fpaths)]), collapse = ", "))
    return(invisible(NULL))
  }
  full <- lapply(fpaths, parse_gsa); names(full) <- TRAITS
  if (any(vapply(full, is.null, logical(1)))) {
    say("SKIP ", coll, " -- unparseable full-data file(s)")
    diag_add("  SKIPPED: unparseable full-data file(s)")
    return(invisible(NULL))
  }
  diag_add("  name column used: ",
           paste(unique(vapply(full, function(v) attr(v, "name_col"), character(1))),
                 collapse = "/"))
  
  jk      <- vector("list", NJK)
  reps_ok <- logical(NJK)
  for (b in seq_len(NJK)) {
    paths <- vapply(TRAITS, function(tr) jk_path(coll, tr, b), character(1))
    if (!all(file.exists(paths))) {
      diag_add("  replicate ", b, " dropped: missing file(s) ",
               paste(basename(paths[!file.exists(paths)]), collapse = ", "))
      next
    }
    vs <- lapply(paths, parse_gsa); names(vs) <- TRAITS
    if (any(vapply(vs, is.null, logical(1)))) {
      diag_add("  replicate ", b, " dropped: unparseable file(s)")
      next
    }
    jk[[b]] <- vs
    reps_ok[b] <- TRUE
  }
  if (sum(reps_ok) < 2L) {
    say("SKIP ", coll, " -- only ", sum(reps_ok), " usable replicate(s)")
    diag_add("  SKIPPED: only ", sum(reps_ok), " usable replicate(s)")
    return(invisible(NULL))
  }
  
  nm_full  <- lapply(full, names)
  nm_jk    <- unlist(lapply(which(reps_ok), function(b) lapply(jk[[b]], names)),
                     recursive = FALSE, use.names = FALSE)
  universe <- Reduce(intersect, c(nm_full, nm_jk))
  universe <- nm_full[[1L]][nm_full[[1L]] %in% universe]   # stable order
  n_uni    <- length(universe)
  
  diag_add("  gene sets in full ", TRAITS[1], " file : ", length(nm_full[[1L]]))
  diag_add("  gene sets kept for this collection : ", n_uni,
           "  (dropped ", length(nm_full[[1L]]) - n_uni, ")")
  diag_add("  usable replicates : ", sum(reps_ok), " of ", NJK)
  rm(nm_full, nm_jk)
  
  if (n_uni < 1L) {
    say("SKIP ", coll, " -- empty universe")
    diag_add("  SKIPPED: empty universe")
    return(invisible(NULL))
  }
  
  beta <- array(NA_real_, dim = c(n_uni, NT, NJK + 1L))
  beta[, , 1L] <- vapply(TRAITS, function(tr) unname(full[[tr]][universe]),
                         numeric(n_uni))
  for (b in which(reps_ok)) {
    beta[, , 1L + b] <- vapply(TRAITS, function(tr) unname(jk[[b]][[tr]][universe]),
                               numeric(n_uni))
    jk[b] <- list(NULL)
  }
  rm(jk, full); invisible(gc(FALSE))
  
  saveRDS(list(collection = coll, sets = universe, traits = TRAITS,
               reps_ok = reps_ok, beta = beta),
          beta_cache(coll), compress = FALSE)
  saveRDS(list(collection = coll, n_sets = n_uni, reps_ok = reps_ok),
          meta_cache(coll))
  
  say(sprintf("%-45s  sets=%6d  usable reps=%3d  %5.1fs", coll, n_uni,
              sum(reps_ok), as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  invisible(NULL)
}

if (STAGE %in% c("all", "parse")) {
  say("stage                : parse")
  say("base                 : ", BASE)
  say("jackknife replicates : ", NJK)
  say("traits (", NT, ")        : ", paste(TRAITS, collapse = ", "))
  say("collections (", length(COLLECTIONS), ")  : ", paste(COLLECTIONS, collapse = ", "))
  say("cache                : ", CACHEDIR)
  
  for (coll in COLLECTIONS) {
    if (RESUME && file.exists(beta_cache(coll)) && file.exists(meta_cache(coll))) {
      say("skip (resume): ", coll); next
    }
    parse_collection(coll)
  }
  
  if (STAGE == "parse") {
    writeLines(DIAG$lines, DIAGFILE)
    say("wrote ", DIAGFILE)
    say("parse stage complete. Run with --stage combine to pool and rank.")
    quit(save = "no", status = 0L)
  }
}


say("stage                : combine")

metas <- list.files(CACHEDIR, pattern = "\\.meta\\.rds$", full.names = TRUE)
if (!length(metas)) stop("No cached collections in ", CACHEDIR,
                         " -- run --stage parse first.", call. = FALSE)
meta <- lapply(metas, readRDS)
colls <- vapply(meta, `[[`, character(1), "collection")
ord   <- order(colls); meta <- meta[ord]; colls <- colls[ord]


reps_ok <- Reduce(`&`, lapply(meta, `[[`, "reps_ok"))
NJK <- length(reps_ok)          # authoritative: taken from the cache, not --n-jk
B <- sum(reps_ok)
if (B < 2L) stop("Only ", B, " replicate(s) usable in every collection.", call. = FALSE)
slices <- c(1L, 1L + which(reps_ok))          # slice 1 = full data
S      <- length(slices)
N      <- sum(vapply(meta, `[[`, numeric(1), "n_sets"))

dropped_reps <- which(!reps_ok)
say("collections pooled   : ", length(colls))
say("gene sets pooled     : ", N)
say("replicates usable    : ", B, " of ", NJK,
    if (length(dropped_reps)) paste0("  (dropped: ",
                                     paste(utils::head(dropped_reps, 10), collapse = ","),
                                     if (length(dropped_reps) > 10) ",..." else "", ")") else "")

read_beta <- function(coll) {
  o <- readRDS(beta_cache(coll))
  if (!identical(o$traits, TRAITS))
    stop("Cached traits for ", coll, " differ from --traits.", call. = FALSE)
  o
}

Ssum <- matrix(0, nrow = NT, ncol = S)
for (cc in colls) {
  o <- read_beta(cc)
  Ssum <- Ssum + colSums(o$beta[, , slices, drop = FALSE], dims = 1L)
  rm(o)
}
mu <- Ssum / N                                   # NT x S

Qsum <- matrix(0, nrow = NT, ncol = S)
for (cc in colls) {
  o <- read_beta(cc)
  d <- sweep(o$beta[, , slices, drop = FALSE], c(2L, 3L), mu, "-")
  Qsum <- Qsum + colSums(d^2, dims = 1L)
  rm(o, d)
}
sdev <- sqrt(Qsum / (N - 1))
if (any(!is.finite(sdev)) || any(sdev == 0))
  stop("Degenerate pooled SD for some trait/replicate.", call. = FALSE)

diag_add("", "=== POOLED NORMALIZATION (full data) ===")
for (i in seq_len(NT))
  diag_add(sprintf("  %-8s mean BETA_STD = %12.6g   sd = %12.6g",
                   TRAITS[i], mu[i, 1L], sdev[i, 1L]))

parts <- vector("list", length(colls))
for (j in seq_along(colls)) {
  o <- read_beta(colls[[j]])
  z <- sweep(sweep(o$beta[, , slices, drop = FALSE], c(2L, 3L), mu, "-"),
             c(2L, 3L), sdev, "/")
  
  comp <- z[, 1L, ]                              # n x S
  for (i in 2:NT) comp <- comp + z[, i, ]
  comp <- comp / NT
  
  c_full <- comp[, 1L]
  cmat   <- comp[, -1L, drop = FALSE]
  cbar   <- rowMeans(cmat)
  se     <- sqrt(((B - 1) / B) * rowSums((cmat - cbar)^2))
  
  parts[[j]] <- data.frame(COLLECTION = colls[[j]], GENE_SET = o$sets,
                           MEAN_Z = c_full, SE = se, stringsAsFactors = FALSE)
  rm(o, z, comp, cmat); invisible(gc(FALSE))
}

all_res <- do.call(rbind, parts); rm(parts)

bad <- !is.finite(all_res$MEAN_Z) | !is.finite(all_res$SE)
if (any(bad)) {
  diag_add("  WARN dropped ", sum(bad), " gene set(s) with non-finite MEAN_Z or SE")
  all_res <- all_res[!bad, , drop = FALSE]
}

dupname <- all_res$GENE_SET[duplicated(all_res$GENE_SET)]
if (length(dupname)) {
  diag_add("", "=== DUPLICATE GENE SET NAMES ACROSS COLLECTIONS ===")
  diag_add("  ", length(unique(dupname)), " name(s) appear in more than one ",
           "collection; each is kept as a separate row (separate MAGMA run) ",
           "and counted once per occurrence in the pooled normalization.")
  diag_add("  e.g. ", paste(utils::head(unique(dupname), 5L), collapse = ", "))
}

if (RANK_WC) {
  all_res <- all_res[order(all_res$COLLECTION, -all_res$MEAN_Z), , drop = FALSE]
  all_res$RANK <- stats::ave(-all_res$MEAN_Z, all_res$COLLECTION,
                             FUN = function(x) seq_along(x))
} else {
  all_res <- all_res[order(-all_res$MEAN_Z), , drop = FALSE]
  all_res$RANK <- seq_len(nrow(all_res))
}

out <- data.frame(
  COLLECTION = all_res$COLLECTION,
  GENE_SET   = all_res$GENE_SET,
  MEAN_Z     = sprintf("%.6f", all_res$MEAN_Z),
  SE         = sprintf("%.6f", all_res$SE),
  RANK       = as.integer(all_res$RANK),
  stringsAsFactors = FALSE
)
utils::write.table(out, OUTFILE, sep = "\t", quote = FALSE,
                   row.names = FALSE, col.names = TRUE)

dfs <- list.files(DIAGDIR, pattern = "\\.txt$", full.names = TRUE)
if (length(dfs))
  DIAG$lines <- c(unlist(lapply(dfs, function(f) readLines(f, warn = FALSE)),
                         use.names = FALSE), DIAG$lines)

absent <- setdiff(COLLECTIONS, colls)
if (length(absent))
  say("WARNING: no cache for ", length(absent), " collection(s): ",
      paste(absent, collapse = ", "))

sp <- split(all_res, all_res$COLLECTION)
lvl <- data.frame(
  COLLECTION = names(sp),
  N          = vapply(sp, nrow, integer(1)),
  MEAN       = vapply(sp, function(d) mean(d$MEAN_Z), numeric(1)),
  SD         = vapply(sp, function(d) stats::sd(d$MEAN_Z), numeric(1)),
  MAX        = vapply(sp, function(d) max(d$MEAN_Z), numeric(1)),
  BEST_RANK  = vapply(sp, function(d) min(d$RANK), numeric(1)),
  stringsAsFactors = FALSE
)
lvl <- lvl[order(-lvl$MEAN), , drop = FALSE]

diag_add("", "=== COLLECTIONS ON THE POOLED SCALE (by mean MEAN_Z) ===")
diag_add(sprintf("  %-42s %7s %9s %8s %9s %10s",
                 "COLLECTION", "N", "MEAN", "SD", "MAX", "BEST_RANK"))
for (r in seq_len(nrow(lvl)))
  diag_add(sprintf("  %-42s %7d %9.4f %8.4f %9.4f %10d",
                   lvl$COLLECTION[r], lvl$N[r], lvl$MEAN[r], lvl$SD[r],
                   lvl$MAX[r], as.integer(lvl$BEST_RANK[r])))

diag_add("", "=== SUMMARY ===")
diag_add("  collections pooled  : ", length(colls))
if (length(absent)) diag_add("  MISSING collections : ", paste(absent, collapse = ", "))
diag_add("  gene sets written   : ", nrow(out))
diag_add("  replicates used (B) : ", B, " of ", NJK)
diag_add("  standardization     : pooled over all ", N, " gene sets in all collections")
diag_add("  ranking             : ",
         if (RANK_WC) "within collection" else "global (all collections pooled)")
diag_add("  finished            : ", format(Sys.time()))
writeLines(DIAG$lines, DIAGFILE)

say("wrote ", OUTFILE, "  (", nrow(out), " gene sets from ", length(colls),
    " collections)")
say("wrote ", DIAGFILE)