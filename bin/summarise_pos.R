require(dplyr)
require(purrr)
require(Rsamtools)

calc.maxpos <- function(id, time, align, posFile, sRNAreads) {
  require(stringr)
  mapInfo <- c("rname", "strand", "pos")
  mapParams <- ScanBamParam(what = mapInfo, tag = c("TC", "TN"))
  bam <- scanBam(file.path(file.home, align), param = mapParams)
  # Now this will ONLY handle files that have tags TC and TN, too!
  map.r <- dplyr::bind_cols(do.call(dplyr::bind_cols, bam[[1]][mapInfo]),
                            do.call(dplyr::bind_cols, bam[[1]]$tag))
  # map.r <- dplyr::bind_cols(bam)
  # Get only + mapping reads. Probably could be done in `scanBam`, too.
  pos.only <-
    map.r %>%
    dplyr::filter(strand == "+")

  r.summary <-
    pos.only %>%
    group_by(rname, pos) %>%
    summarise(count = n()) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname))

  tc.summary <-
    pos.only %>%
    group_by(rname, pos) %>%
    filter(!is.na(TC)) %>%
    summarise(tcReads = n()) %>% ungroup() %>%
    mutate(flybase_id = as.character(rname)) %>%
    select(-rname)

  r.sum.pos <-
    r.summary %>%
    left_join(tc.summary, by = c("flybase_id", "pos")) %>%
    left_join(mir.anno, by = "flybase_id") %>%
    select(-rname, -loop)

  r.sum.arms <-
    r.sum.pos %>%
    dplyr::filter((pos >= `5p` - 5 & pos <= `5p` + 5) | (pos >= `3p` - 5 & pos <= `3p` + 5)) %>%
    mutate(arm.name = ifelse(pos >= `5p` - 5 & pos <= `5p` + 5,
                             paste0(str_sub(mir_name, 5, -1), "-5p"),
                             paste0(str_sub(mir_name, 5, -1), "-3p")))

  countName <- paste("totalReads", id, sep = ".")
  tcName <- paste("tcReads", id, sep = ".")

  r.sum.max.pos <-
    r.sum.arms %>%
    group_by(arm.name) %>%
    top_n(n = 5, wt = count) %>% ungroup() %>%
    mutate(count = count / sRNAreads * 1000000,
           tcReads = tcReads / sRNAreads * 1000000) %>%
    rename_(.dots = setNames(c("count", "tcReads"), c(countName, tcName)))

  return(r.sum.max.pos)
}
