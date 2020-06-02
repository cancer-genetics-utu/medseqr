convertId2 <-
		function (id, species = "Human")
{
	if (species == "Human") {
		ensg2eg.env <- org.Hs.egENSEMBL2EG
		sym.env <- org.Hs.egSYMBOL
		sym2eg.env <- org.Hs.egSYMBOL2EG
		ensg.env <- org.Hs.egENSEMBL
		ensg <- "ENSG"
	}
	if (species == "Mouse") {
		ensg2eg.env <- org.Mm.egENSEMBL2EG
		sym.env <- org.Mm.egSYMBOL
		sym2eg.env <- org.Mm.egSYMBOL2EG
		ensg.env <- org.Mm.egENSEMBL
		ensg <- "ENSMU"
	}
	if (length(id) == 1) {
		if (length(grep(ensg, id)) > 0) {
			if (AnnotationDbi::exists(id, envir = ensg2eg.env)) {
				entrez <- get(id, envir = ensg2eg.env)
				if (length(entrez) > 1) {
					sym <- NA_character_
				}
				else {
					if (AnnotationDbi::exists(entrez, envir = sym.env)) {
						sym <- paste(get(entrez, envir = sym.env),
										collapse = " /// ")
					} else {
						sym <- NA_character_
					}
				}
			}
			else {
				sym <- NA_character_
			}
		}
		else {
			if (AnnotationDbi::exists(id, envir = sym2eg.env)) {
				entrez <- get(id, envir = sym2eg.env)
				if (length(entrez) > 1) {
					sym <- NA_character_
				}
				else {
					if (AnnotationDbi::exists(entrez, envir = ensg.env)) {
						sym <- paste(get(entrez, envir = ensg.env),
										collapse = " /// ")
					} else {
						sym <- NA_character_
					}
				}
			}
			else {
				sym <- NA_character_
			}
		}
		names(sym) <- id
		return(sym)
	}
	else {
		if (length(grep(ensg, id[1])) > 0) {
			entrez <- mget(id, envir = ensg2eg.env, ifnotfound = NA)
			entrez <- sapply(entrez, function(x) {
						if (length(x) > 1 || is.na(x)) {
							"---"
						}
						else {
							x
						}
					})
			hugo <- mget(entrez, envir = sym.env, ifnotfound = NA)
			hugo <- sapply(hugo, function(x) {
						if (length(x) > 1) {
							paste(x, collapse = " /// ")
						}
						else {
							x
						}
					})
			names(hugo) <- id
			return(hugo)
		}
		else {
			entrez <- mget(id, envir = sym2eg.env, ifnotfound = NA)
			entrez <- sapply(entrez, function(x) {
						if (length(x) > 1 || is.na(x)) {
							"---"
						}
						else {
							x
						}
					})
			ensg <- mget(entrez, envir = ensg.env, ifnotfound = NA)
			ensg <- sapply(ensg, function(x) {
						if (length(x) > 1) {
							paste(x, collapse = " /// ")
						}
						else {
							x
						}
					})
			names(ensg) <- id
			return(ensg)
		}
	}
}

#' @export
convert.bm <-
		function(dat, id="ID", biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride"),
				host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"),
				rm.dups=FALSE)
{
	if (id=="row.names") {
		values <- rownames(dat)
	} else {
		values <- dat[[id]]
	}
	biom.ids <- get.bm(values, biom.data.set, biom.mart, host, biom.filter, biom.attributes)
	gene.lab <- merge(biom.ids, dat, by.x=biom.filter, by.y=id, all.y=TRUE, all.x=FALSE, sort=TRUE)
	if (rm.dups) {
		cat("  Removing", length(which(duplicated(gene.lab[[biom.filter]]))), "duplicated row(s)...\n")
		gene.lab <- gene.lab[!duplicated(gene.lab[[biom.filter]]), ]
	}
	if (any(gene.lab$hgnc_symbol=="") || any(is.na(gene.lab$hgnc_symbol))) {
		cat("  Replacing", length(which(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))), "missing Gene Symbols by Ensembl IDs...\n")
		gene.lab$hgnc_symbol[gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol)] <- gene.lab[[biom.filter]][gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol)]
	}
	return(gene.lab)
}

# fetch additional information from biomart
#' @export
get.bm <-
		function(values, biom.data.set="hsapiens_gene_ensembl", biom.mart=c("ensembl", "snp", "funcgen", "vega", "pride"),
				host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol","description"))
{
	marts <- listMarts(host=host)[["biomart"]]
	marts1 <- unlist(lapply(strsplit(tolower(marts), "_"), function(x) x[length(x)]))
	biom <- match.arg(biom.mart)
	biom <- marts[grep(biom, marts1)]
	cat("Using BioMart:", biom, "\n")
	mart <- useDataset(dataset=biom.data.set, mart=useMart(biomart=biom, host=host))
	cat("  Getting more information:", biom.attributes[biom.attributes!=biom.filter], "...\n")
	getBM(attributes=biom.attributes, filters=biom.filter, values=as.character(values), mart=mart)
}

#' @export
todisp2 <-
		function(ensg, lab=NULL, biomart=TRUE)
{
	if (biomart) {
		sym <- get.bm(ensg, biom.data.set="hsapiens_gene_ensembl", biom.mart="ensembl", host="www.ensembl.org", biom.filter="ensembl_gene_id", biom.attributes=c("ensembl_gene_id","hgnc_symbol"))
	} else if(!is.null(lab)) {
		cat("  Using input data frame for ID conversion...\n")
		sym <- data.frame(ensembl_gene_id=rownames(lab), hgnc_symbol=lab[, 1], stringsAsFactors=FALSE)
	} else {
		cat("  Using 'AnnotationDbi framework for ID conversion...\n")
		sym <- convertId2(ensg)
		return(ifelse(is.na(sym), ensg, sym))
	}
	gene.lab <- merge(data.frame(ensembl_gene_id=ensg, stringsAsFactors=FALSE), sym, by="ensembl_gene_id", sort=FALSE)
	if (any(gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol))) {
		cat("    Replacing", length(which(gene.lab$hgnc_symbol=="")), "missing Gene Symbols by Ensembl IDs...\n")
		replace <- gene.lab$hgnc_symbol=="" | is.na(gene.lab$hgnc_symbol)
		gene.lab$hgnc_symbol[replace] <- gene.lab$ensembl_gene_id[replace]
	}
	gene.lab$hgnc_symbol
}
