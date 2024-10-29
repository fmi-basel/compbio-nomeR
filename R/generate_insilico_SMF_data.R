#' Generate in-silico single-molecule footprinting dataset
#'
#' @param region_len length of the region
#' @param n_reads number of molecule to generate
#' @param footprint_models \code{list} describing footprints which should be generated. Each footprint within footprint_models is also a \code{list} that must contain
#' \code{PROTECT_PROB} - vector of emission probabilities with same length and footprint length;
#' \code{COVER_PROB} - \code{numeric} that describes footprint coverage
#' \code{NAME} - \code{character} with the name of the footprint
#' @param ftp_emission_posprob \code{matrix} describing positional probabilities for each footprint. Number of rows must be equal to number of footprints and number of columns must be equal to \code{region_len}
#' @param bgprotectprob \code{numeric} emission probability for background
#' @param infposdens if scalar between 0 and 1 treated as percentage of informative positions.
#' If a vector of integers treated as predefined informative positions.
#'
#' @return \code{list} that contains
#' "DATA" - \code{matrix} with in-silico data;
#' "TRUE_CONF" - true positions for each footprint;
#' "ftp_pos_prob" -normalized probabilities of placing each footprint.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ## length of a chromosome
#' chr_len <- 1000
#'
#' ftp_cover_c2 <- c(0.05,0.5)
#' ftp_len_c2 <- c(50,150)
#' case2_pos_mat <- rbind(matrix(c(rep(0,5150),1,rep(0,4849)),nrow=1,ncol=chr_len),
#' 											 matrix(1/chr_len,nrow=1,ncol=chr_len))
#'
#' ## sequencing depths
#' max_nreads <- 30
#' ## footprint model
#' case2_ftp_models <- list("ftp50" = list("PROTECT_PROB" = rep(1,ftp_len_c2[1]),
#' 																				"COVER_PRIOR" = ftp_cover_c2[1],
#' 																				"NAME" = "ftp50"),
#' 												 "ftp150"=list("PROTECT_PROB" = rep(1,ftp_len_c2[2]),
#' 												 							"COVER_PRIOR" = ftp_cover_c2[2],
#' 												 							"NAME" = "ftp150"))
#'
#' smf_data <- generate_insilico_SMF_data(region_len = chr_len, # length of the region
#' 																			n_reads = max_nreads, # number of reads
#' 																			footprint_models = case2_ftp_models,
#' 																			ftp_emission_posprob=case2_pos_mat, ## matrix with nrow=number ftps ncol=ampllen for each footprint to describe
#' 																			## probabilities of emission for each amplicon position
#' 																			bgprotectprob = 0,
#' 																			infposdens = 1 # if scalar between 0 and 1 treated as percentage of informative positions
#' 																			# if vector of integers treated as predefined informative positions
#' 																			)
#'
#'}
#'
#'
#'
generate_insilico_SMF_data <- function(region_len, # length of the amplicon
																			 n_reads, # number of reads
																			 footprint_models,
																			 ftp_emission_posprob=NULL, ## matrix with nrow=number ftps ncol=ampllen for each footprint to describe
																			 ## probabilities of emission for each amplicon position

																			 bgprotectprob,
																			 infposdens = 1
){


	extend_ampl = FALSE
	#### CHECK ARGUMENTS ####
	checkmate::assert_integerish(region_len,len = 1,any.missing = F)
	checkmate::assert_integerish(n_reads,len = 1,any.missing = F)
	if(checkmate::testNumber(infposdens,lower = 0.0,upper = 1.0, null.ok = F)){
		checkmate::assert_number(infposdens,lower = 0.0,upper = 1.0, null.ok = F)
	} else {
		checkmate::assert_integerish(infposdens,lower=1,upper=region_len,any.missing = F,min.len = 1,unique = T,null.ok = F)
	}
	checkmate::assert_list(footprint_models,types = "list",any.missing = F,min.len = 1)
	lapply(footprint_models,
				 function(x){
				 	checkmate::assert_numeric(x$PROTECT_PROB,lower = 0.0,upper = 1.0,finite = T,null.ok = F,any.missing = F,min.len = 2)
				 	checkmate::assert_number(x$COVER_PRIOR,lower = 0.0,upper = 1.0,finite = T,null.ok = F)
				 	checkmate::assert_character(x$NAME,min.chars = 1,any.missing = F,len=1)
				 })

	checkmate::assert_matrix(ftp_emission_posprob,
													 mode = "numeric",
													 any.missing = F,
													 nrows = length(footprint_models),
													 ncol=region_len,
													 null.ok = T)

	checkmate::assert_number(bgprotectprob,lower = 0.0,upper = 1.0,null.ok = F,finite = T)

	#### GET INFORMATIVE POSITIONS ####

	if(length(infposdens) == 1 & all(infposdens > 0) & all(infposdens <= 1)){
		inf_pos <- sort(sample(1:region_len,size=floor(infposdens * region_len),
													 replace=F)) # replace = F makes sure that positions will not repeat

	} else if(is.vector(infposdens,"integer") | is.vector(infposdens,"numeric")){
		infposdens <- as.integer(infposdens)
		if(any(!infposdens %in% 1:region_len)){
			warnings("infposdens contains positions outsite amplicon length. Removing such positions")
		}
		inf_pos <- infposdens[infposdens %in% 1:region_len]
	} else {
		stop(paste0("infposdens must be a scalar between 0 and 1 or a vector of integers or numerics."))
	}

	#browser()
	#### PROCESSING FOOTPRINT MODELS #####
	## check whether any of the footprint has length 1 or name "BG"
	if(any(sapply(footprint_models,function(x){
		x$NAME
	}) == "BG") | any(sapply(footprint_models,function(x){
		length(x$PROTECT_PROB)
	}) == 1)){
		stop("Footprint with name BG or length 1 are reserved for background. Remove them from footprint models.")
	}


	## adding BG
	bgcoverprior <- 1 - sum(sapply(footprint_models,function(x){
		x$COVER_PRIOR
	},simplify = T,USE.NAMES = F))

	footprint_models <- c(list(list("NAME" = "BG",
																	"COVER_PRIOR" = bgcoverprior,
																	"PROTECT_PROB" = bgprotectprob)),
												footprint_models)

	#browser()
	model_names <- sapply(footprint_models,function(x){
		x$NAME
	})

	if(any(duplicated(model_names))){
		stop("Footprint names must not be duplicated!")
	}

	names(footprint_models) <- model_names

	model_cover_priors <- sapply(footprint_models,function(x){
		x$COVER_PRIOR
	},simplify = T,USE.NAMES = F)

	model_lengths <- sapply(footprint_models,function(x){
		length(x$PROTECT_PROB)
	},simplify = T,USE.NAMES = F)


	model_start_priors <- .cover_prior2start_prior(model_cover_priors,
																								 model_lengths)

	#browser()

	## decide whether start of footprints will be outside amplicon
	if(extend_ampl){
		max.model.len <- 2*max(model_lengths)
	} else{
		max.model.len <- 1
	}

	reglen_ext <- region_len + max.model.len - 1
	#browser()

	### assign footprint start probabilities for each position of the amplicon
	if(is.null(ftp_emission_posprob)){
		ftp_emission_posprob <- matrix(data = 1/reglen_ext,
																	 nrow=length(model_names),
																	 ncol = reglen_ext)
	} else {
		## add first row for background
		ftp_emission_posprob <- rbind(rep(1/region_len,region_len),
																	ftp_emission_posprob)

		## if start of ftps outside amplicon boundaries allowed, add additional probs to the left and renormalize
		if(extend_ampl){
			ext_mat <- matrix(data = ftp_emission_posprob[,1],
												nrow = nrow(ftp_emission_posprob),
												ncol = max.model.len - 1,byrow=FALSE)
			ftp_emission_posprob <- cbind(ext_mat,
																		ftp_emission_posprob)
		} else {
			## if extention of amplicon is not allowed, set positional probabilities of footprints which do not fit at some positions to 0
			ftp_emission_posprob <- do.call(rbind,lapply(1:nrow(ftp_emission_posprob),
																									 function(ftp){
																									 	mod_len <- model_lengths[ftp]
																									 	prvec <- ftp_emission_posprob[ftp,]
																									 	zero_start_pos <- reglen_ext - mod_len+2
																									 	if(zero_start_pos <= reglen_ext){
																									 		prvec[(reglen_ext-mod_len+1):reglen_ext] <- 0
																									 	}
																									 	return(prvec)
																									 }))
		}
	}




	ftp_emission_posprob <- ftp_emission_posprob / rowSums(ftp_emission_posprob)

	ftp_emission_posprob <- apply(ftp_emission_posprob,MARGIN = 2,
																function(x){
																	x*model_start_priors/sum(x*model_start_priors)
																})
	#browser()
	#require(parallel)
	SMF_data <- lapply(1:n_reads,function(i){

		dat <- matrix(NA,nrow = 1,ncol = reglen_ext)
		true_conf_data <- data.frame(seq = i,
																 pos = 1:reglen_ext - max.model.len + 1,
																 model = NA_character_,
																 stringsAsFactors = F)
		# true_conf_pos <- vector(mode="integer")
		# true_conf_model <- vector(mode = "character")
		pos <- 1
		while(pos <= reglen_ext){

			## choose which model to emit
			emit_model <- sample(x = model_names,
													 size = 1,
													 prob = ftp_emission_posprob[,pos])

			#true_conf_pos <- c(true_conf_pos,pos - max.model.len + 1)
			# true_conf_pos <- c(true_conf_pos,pos)
			# true_conf_model <- c(true_conf_model,emit_model)
			true_conf_data$model[pos] <- emit_model

			model_len <- model_lengths[emit_model]
			model_prot_prob <- footprint_models[[emit_model]][["PROTECT_PROB"]]
			## draw random data using binomial distribution

			model_synth_data <- rbinom(n = model_len, size = 1, prob = model_prot_prob)
			fill.idx <- pos:(min(pos + model_len - 1,reglen_ext))
			dat[1,fill.idx] <- model_synth_data[1:length(fill.idx)]
			pos <- pos + model_len
			names(pos) <- NULL

		}



		dat.out <- list("DATA" = dat,
										"TRUE_CONF" = subset(true_conf_data,!is.na(model)))
		return(dat.out)
	})

	SMF_sim_data <- do.call(rbind,lapply(SMF_data,function(x){
		x[["DATA"]]
	}))
	true_conf <- do.call(rbind,lapply(SMF_data,function(x){
		x[["TRUE_CONF"]]
	}))
	## remove additional flanks
	SMF_sim_data <- SMF_sim_data[,max.model.len:reglen_ext,drop=F]
	## set non GpC positions to NA
	SMF_sim_data[,!(1:ncol(SMF_sim_data) %in% inf_pos)] <- NA

	## return without additional flanking columns
	list("DATA" = SMF_sim_data,
			 "TRUE_CONF" = true_conf,
			 "ftp_pos_prob" = ftp_emission_posprob)
}


