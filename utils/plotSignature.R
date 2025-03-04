toPyr <- function(nctds) {
  nctds.pyr <- nctds
  nctds.pyr[nctds=='A'] <- 'T'
  nctds.pyr[nctds=='T'] <- 'A'
  nctds.pyr[nctds=='G'] <- 'C'
  nctds.pyr[nctds=='C'] <- 'G'
  as.character(nctds.pyr)
}

get_context = function(subs.df) {
  subs.df$bb <- substr(subs.df$triplets,1,1)
  subs.df$ba <- substr(subs.df$triplets,3,3)
  subs.df$isPyr <- subs.df$Ref %in% c('C', 'T')
  subs.df$context96 <- ''
  subs.df$context96 <- paste0(subs.df$bb,'[', subs.df$Ref,'>',subs.df$Alt,']', subs.df$ba )
  subs.df$context96[!subs.df$isPyr] <- paste0(toPyr(subs.df$ba[!subs.df$isPyr]),
                                              '[', toPyr(subs.df$Ref[!subs.df$isPyr]),'>',
                                              toPyr(subs.df$Alt[!subs.df$isPyr]),']', 
                                              toPyr(subs.df$bb[!subs.df$isPyr]))
  
  subs.df = subs.df %>% dplyr::select(-bb, -ba, -isPyr)
  return(subs.df)
}

get_context_table = function(subs.df, mut.order) {
  preBase <- rep(c(rep('A', 4), rep('C', 4), rep('G', 4), rep('T', 4)),6)
  refBase <- c(rep('C', 48 ), rep('T', 48))
  altBase <- c(rep('A', 16 ), rep('G', 16), rep('T', 16),rep('A', 16 ), rep('C', 16), rep('G', 16))
  postBase <- rep(c('A', 'C', 'G', 'T'), 96/4)
  mut.order <- paste0(preBase, '[',refBase, '>', altBase, ']', postBase )
  
  b <- table(subs.df$context96)
  # fill 0
  b[mut.order[!(mut.order %in% names(b))]] = 0
  b <- b[as.character(mut.order)]	
  b[is.na(b)] <- 0	
  
  b = b %>% as.matrix %>% as.data.frame %>%
    tibble::rownames_to_column('context96') %>%
    dplyr::rename(n = V1) %>%
    mutate(substitution = substring(context96, 3, 5)) %>%
    mutate(context96 = factor(context96, mut.order)) %>%
    filter(!is.na(context96))
  
  return(b)
}

signature_colors = c(
  'C>A' = 'cyan3',
  'C>G' = 'black',
  'C>T' = 'red',
  'T>A' = 'grey',
  'T>C' = 'darkolivegreen3',
  'T>G' = 'hotpink'
)

plot_signature = function(context_table, title = '') {
  
  p = ggplot(context_table) +
    geom_col(aes(x = context96, y = n, fill = substitution), color = 'black', size = 0.3) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    theme(plot.title = element_text(hjust = 0.5, size = 20)) +
    scale_fill_manual(values = signature_colors) +
    ggtitle(paste(title, '\n', 'n = ', sum(context_table$n)))
  
  return(p)
}