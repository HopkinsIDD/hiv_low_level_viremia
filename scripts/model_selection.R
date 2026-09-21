suppressMessages(library(loo))
suppressMessages(library(cmdstanr))
suppressMessages(library(tidyverse))


fit_paths = c(
	'output/ve_constrained',
	'output/ve_constrained_estC',
	'output/ve_nosaturation',
	'output/sigmoid_constrained')
labels = c('main_analysis',
	'estC',
	'no_saturation',
	'sigmoid')



loos = list()
for (idx in seq(1, length(fit_paths))){
	fit_path = fit_paths[idx]
	# read in periods
	periods = read_tsv('data/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv', show_col_types=FALSE) %>%
		mutate(period_idx = as.character(seq(1,n())))
	#copies = do.call(
	#	rbind,
	#	lapply(
	#		strsplit(
	#			gsub("\\.", '', gsub('\\]', "", gsub('\\[', "", periods$donor_copies))), 
	#			"\\s+"), 
	#		function(x) {
	#	  		length(x) <- 2
	#	  		return(x)}))
	#periods$donor_copies1 = as.numeric(copies[,1])
	#periods$donor_copies2 = as.numeric(copies[,2])
	#periods = periods %>% mutate(min_copies = pmin(donor_copies1, donor_copies2, na.rm=TRUE))
	fit_draws = read_tsv(paste0(fit_path, '/log_prob_no_transmission.tsv'))
	# get no likelihood of no intra-couple transmission
	log_likelihood = fit_draws[colnames(fit_draws)[grepl('log_prob_no_transmission\\[', colnames(fit_draws))]] %>%
		mutate(cit = seq(1, n()))%>%
		pivot_longer(-c(cit)) %>%
		mutate(
			period_idx = gsub('\\]', '', str_split(name, '\\[', simplify=TRUE)[,2])) %>%
		left_join(periods, by='period_idx') %>%
		mutate(
			log_likelihood = if_else(
				conversion == FALSE,
				value,
				log(1 - exp(value)))) %>%
		pivot_wider(id_cols=cit, names_from=period_idx, values_from=log_likelihood) %>%
		select(-cit)
	loo_fit = loo(as.matrix(log_likelihood), save_psis = TRUE)
	loos[[labels[idx]]] = loo_fit
	#loo_fit2 = loo(as.matrix(log_likelihood2), save_psis = TRUE)
	#loos2[[labels[idx]]] = loo_fit2
}


comp = loo_compare(loos)
print(comp, digits=3, p_worse=TRUE)

write_tsv(
	as_tibble(comp),
	'output/transmission_model_selection.tsv')
