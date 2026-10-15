feats_df <- function(...) {
  rows <- list(...)
  data.frame(path = vapply(rows, `[[`, 1, 1), scaffold = vapply(rows, `[[`, 1, 2),
    element = vapply(rows, `[[`, "", 3), stringsAsFactors = FALSE)
}

test_that("contigs sharing an anchor form one set; stray contigs are dropped", {
  f <- feats_df(list(1, 1, "18S"), list(1, 1, "28S"), list(2, 1, "18S"), list(2, 2, "28S"),
    list(3, 1, "ITS1"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(s[[1]]$units, data.frame(path = c(1L, 2L, 2L), scaffold = c(1L, 1L, 2L)), ignore_attr = TRUE)
  expect_equal(s[[1]]$label, "18S-28S (3 contigs)")
})

test_that("three single-contig paths are one set of three", {
  f <- feats_df(list(1, 1, "18S"), list(2, 1, "18S"), list(3, 1, "18S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 1L)
  expect_equal(nrow(s[[1]]$units), 3L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
})

test_that("disjoint sets come out largest first and singletons vanish", {
  f <- feats_df(list(1, 1, "28S"), list(2, 1, "28S"), list(1, 2, "18S"), list(2, 2, "18S"),
    list(3, 2, "18S"), list(4, 1, "5.8S"))
  s <- RiboPilot:::paths_homology_sets(f)
  expect_length(s, 2L)
  expect_equal(s[[1]]$label, "18S (3 contigs)")
  expect_equal(s[[2]]$label, "28S (2 contigs)")
  expect_equal(s[[2]]$units$path, c(1L, 2L))
})

test_that("no anchors or empty input gives an empty list", {
  expect_equal(RiboPilot:::paths_homology_sets(feats_df(list(1, 1, "ITS1"), list(2, 1, "ITS2"))), list())
  expect_equal(RiboPilot:::paths_homology_sets(data.frame(path = integer(), scaffold = integer(),
    element = character())), list())
})

test_that("a link function vetoes a shared anchor, but a third contig still joins both", {
  f <- feats_df(list(1, 1, "28S"), list(2, 1, "28S"), list(3, 1, "5.8S"))
  veto <- function(ua, ub, el) !(el == "28S" && setequal(c(ua, ub), c("1 1", "2 1")))
  expect_equal(RiboPilot:::paths_homology_sets(f, link = veto), list())
  f <- feats_df(list(1, 1, "28S"), list(2, 1, "28S"), list(1, 1, "18S"), list(2, 1, "18S"),
    list(3, 1, "18S"))
  s <- RiboPilot:::paths_homology_sets(f, link = function(ua, ub, el)
    veto(ua, ub, el) && !(el == "18S" && setequal(c(ua, ub), c("1 1", "2 1"))))
  expect_length(s, 1L)
  expect_equal(s[[1]]$units$path, c(1L, 2L, 3L))
})

test_that("the overlap link joins copies that align and rejects ones that do not", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  set.seed(7)
  a <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
  b <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ? WHERE path = 1", params = list(a))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ? WHERE path = 2 AND scaffold = 1", params = list(b))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'28S',1,300,'+'), ('S1',1,2,'28S',1,300,'+'), ('S1',2,1,'28S',1,300,'+'),
    ('S1',2,2,'28S',1,30,'+')")
  link <- RiboPilot:::paths_overlap_link(con, "S1")
  DBI::dbDisconnect(con)
  expect_true(link("1 1", "1 2", "28S"))
  expect_true(link("1 2", "1 1", "28S"))
  expect_false(link("1 1", "2 1", "28S"))
  expect_false(link("1 1", "2 2", "28S"))   # piece shorter than REF_ALN_MIN_BP
  expect_false(link("1 1", "1 2", "18S"))   # no such element
})

aln_fixture <- function() {
  Biostrings::DNAStringSet(c(`p1.s1` = "ACGTACGTAC", `p2.s1` = "ACG--CGTAC"))
}
no_feats <- function() data.frame(unit = character(), element = character(), start = numeric(),
  stop = numeric(), strand = character(), partial_start = numeric(), partial_stop = numeric())

test_that("the payload maps features from contig to alignment columns across gaps", {
  f <- data.frame(unit = c("p1.s1", "p2.s1", "p2.s1"), element = c("18S", "18S", "ITS1"),
    start = c(2, 2, NA), stop = c(5, 5, 7), strand = "+", partial_start = c(0, 1, 0),
    partial_stop = c(0, 0, 0), stringsAsFactors = FALSE)
  p <- RiboPilot:::paths_payload(aln_fixture(), f, c(p1.s1 = "p1.s1 keep", p2.s1 = "p2.s1 drop"))
  expect_equal(p$len, 10L)
  expect_length(p$rows, 2L)
  expect_equal(p$rows[[2]]$label, "p2.s1 drop")
  expect_equal(p$rows[[2]]$seq, "ACG--CGTAC")
  r1 <- p$rows[[1]]$features[[1]]; r2 <- p$rows[[2]]$features[[1]]
  expect_equal(c(r1$pos1, r1$pos2), c(2L, 5L))
  expect_equal(c(r2$pos1, r2$pos2), c(2L, 7L))   # bases 2..5 of p2.s1 sit at columns 2,3,6,7
  expect_true(r2$partial5); expect_false(r1$partial5)
  expect_equal(r1$color, unname(RiboPilot:::RRNA_ELEMENT_COLORS["18S"]))
  expect_length(p$rows[[2]]$features, 1L)        # NA start is skipped
})

test_that("the difference count includes gap columns and is zero for identical rows", {
  lab <- c(a = "a", b = "b", p1.s1 = "a", p2.s1 = "b")
  rows_of <- function(x) RiboPilot:::paths_payload(x, no_feats(), lab)$rows
  expect_equal(RiboPilot:::paths_diff_count(rows_of(aln_fixture())), 2L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGT", b = "ACGT")))), 0L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGTAC", b = "--GTAC")))), 0L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGTAC", b = "AC-TAC")))), 1L)
  expect_equal(RiboPilot:::paths_diff_count(rows_of(Biostrings::DNAStringSet(c(a = "ACGTACGT", b = "--GTTC--")))), 1L)
})

test_that("paths_align aligns two near-identical sequences to equal width", {
  skip_if_not_installed("DECIPHER")
  a <- strrep("ACGTTGCAAC", 8)
  b <- paste0(substr(a, 1, 40), substr(a, 44, 80))
  m <- RiboPilot:::paths_align(c(`p1.s1` = a, `p2.s1` = b))
  expect_s4_class(m, "DNAStringSet")
  expect_equal(names(m), c("p1.s1", "p2.s1"))
  expect_equal(length(unique(Biostrings::width(m))), 1L)
  expect_null(RiboPilot:::paths_align(c(`p1.s1` = a)))
})

test_that("paths_send aligns the set's contigs and sends rows to the paths canvas", {
  skip_if_not_installed("withr"); skip_if_not_installed("DECIPHER")
  dir <- mk_app_project(); withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  a <- strrep("ACGTTGCAAC", 8); b <- paste0(substr(a, 1, 40), substr(a, 44, 80))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, triage = 'keep' WHERE path = 1 AND scaffold = 1", params = list(a))
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?, triage = 'drop' WHERE path = 2 AND scaffold = 1", params = list(b))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,30,'+'), ('S1',2,1,'18S',1,30,'+')")
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  set <- list(label = "18S (2 contigs)", units = data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  p <- RiboPilot:::paths_send(session, "S1", set)
  expect_equal(sent$id, "seqv_paths-canvas")
  expect_equal(sent$unit, "S1")
  expect_equal(vapply(sent$rows, `[[`, "", "label"), c("p1.s1 keep", "p2.s1 drop"))
  expect_equal(sent$rows[[1]]$features[[1]]$gene, "18S")
  expect_identical(p$rows, sent$rows)
})

test_that("paths_labels adds triage to the key and leaves it bare when unset", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET triage = 'keep' WHERE path = 1 AND scaffold = 1")
  DBI::dbExecute(con, "UPDATE assemblies SET triage = NULL WHERE path = 2 AND scaffold = 1")
  lab <- RiboPilot:::paths_labels(con, "S1", data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  DBI::dbDisconnect(con)
  expect_equal(lab, c(p1.s1 = "p1.s1 keep", p2.s1 = "p2.s1"))
})

test_that("paths_send marks a feature that touches the contig end as partial", {
  skip_if_not_installed("withr"); skip_if_not_installed("DECIPHER")
  dir <- mk_app_project(); withr::local_dir(dir)
  set.seed(7)
  a <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ? WHERE scaffold = 1 AND path IN (1, 2)", params = list(a))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'28S',101,300,'+'), ('S1',2,1,'28S',101,300,'+')")
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  set <- list(label = "28S (2 contigs)", units = data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  RiboPilot:::paths_send(session, "S1", set)
  f <- sent$rows[[1]]$features[[1]]
  expect_true(f$partial3)
  expect_false(f$partial5)
})

test_that("paths_send returns NULL when fewer than two contigs have a sequence", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  con <- RiboPilot:::app_con()
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = NULL WHERE path = 2")
  DBI::dbDisconnect(con)
  sent <- NULL
  session <- list(sendCustomMessage = function(type, msg) sent <<- msg)
  set <- list(label = "x", units = data.frame(path = c(1L, 2L), scaffold = c(1L, 1L)))
  expect_null(RiboPilot:::paths_send(session, "S1", set))
  expect_null(sent)
})

test_that("the viewer's align mode has difference buttons and no coverage boxes", {
  html <- paste(as.character(RiboPilot:::seqview_ui("seqv_paths", NULL, mode = "align")), collapse = " ")
  expect_match(html, 'data-rpseq="prev_diff"', fixed = TRUE)
  expect_match(html, 'data-rpseq="next_diff"', fixed = TRUE)
  expect_match(html, "seqv_paths-canvas", fixed = TRUE)
  expect_match(html, "seqv_paths-show_nt", fixed = TRUE)
  expect_false(grepl("show_cov", html, fixed = TRUE))
  expect_false(grepl('data-rpseq="fit"', html, fixed = TRUE))
  expect_false(grepl("<details", html, fixed = TRUE))
  unit <- paste(as.character(RiboPilot:::seqview_ui("seqv_asm", "Cov")), collapse = " ")
  expect_match(unit, "show_cov", fixed = TRUE)
  expect_false(grepl("prev_diff", unit, fixed = TRUE))
})

test_that("the Assemble window shows the paths section only for a multi-path sample with a set", {
  skip_if_not_installed("withr")
  dir <- mk_app_project(); withr::local_dir(dir)
  html <- function() paste(as.character(RiboPilot:::assemble_modal("S1")), collapse = " ")
  expect_false(grepl("seqv_paths-canvas", html(), fixed = TRUE))   # no annotations yet
  con <- RiboPilot:::app_con()
  set.seed(3)
  DBI::dbExecute(con, "UPDATE assemblies SET sequence = ?",
    params = list(paste(sample(c("A", "C", "G", "T"), 200, TRUE), collapse = "")))
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,1,'18S',1,120,'+'), ('S1',2,1,'18S',1,120,'+')")
  h <- html()
  expect_match(h, "seqv_paths-canvas", fixed = TRUE)
  expect_match(h, "Compare assembly paths", fixed = TRUE)
  expect_false(grepl("paths_set", h, fixed = TRUE))                 # one set: no dropdown
  DBI::dbExecute(con, "INSERT INTO annotations (ID, path, scaffold, element, start, stop, strand)
    VALUES ('S1',1,2,'28S',121,200,'+'), ('S1',2,2,'28S',121,200,'+')")
  expect_match(html(), 'id="paths_set"', fixed = TRUE)
  DBI::dbExecute(con, "DELETE FROM assemblies WHERE path = 2")
  DBI::dbExecute(con, "DELETE FROM annotations WHERE path = 2")
  DBI::dbDisconnect(con)
  expect_false(grepl("seqv_paths-canvas", html(), fixed = TRUE))   # one path left
})

test_that("the summary line counts rows, columns, and differences", {
  rows <- RiboPilot:::paths_payload(aln_fixture(), no_feats(), c(p1.s1 = "a", p2.s1 = "b"))
  expect_equal(RiboPilot:::paths_summary_text(rows), "2 contigs, 10 columns, 2 differ")
  expect_equal(RiboPilot:::paths_summary_text(NULL), "Could not align these contigs.")
  big <- list(len = 7812L, rows = rows$rows)
  expect_match(RiboPilot:::paths_summary_text(big), "7,812 columns", fixed = TRUE)
})

test_that("the align-mode viewer box is height-capped and scrolls", {
  html <- paste(as.character(RiboPilot:::seqview_ui("seqv_paths", NULL, mode = "align")), collapse = " ")
  expect_match(html, 'class="rp-seqview rp-seqview-align"', fixed = TRUE)
  unit <- paste(as.character(RiboPilot:::seqview_ui("seqv_asm", "Cov")), collapse = " ")
  expect_false(grepl("rp-seqview-align", unit, fixed = TRUE))
  expect_equal(RiboPilot:::PATHS_WARN_N, 20L)
})
