
# Script written by AM Senior @ the University of Sydney to generate phylogenetic covariance matrices for various meta-analytic datasets on diet mixing

# Clean up
rm(list=ls())

# Load packages
library(rotl)
library(ape)
library(stringr)

httr::set_config(httr::config(ssl_verifypeer = FALSE))

# Load the data
data<-read.csv("data/full_data_phil_trans.csv")

# Create a copy of the dataset with some species recoded as some taxa as only their sister species can be found on iOTL (note these were found by trial and error) - sister species found by checking publiched phylogenies
# These can be recoded on the tree after
pseudo_data<-data
missing_taxa<-c("Nabis_roseipennis", "Balaustium_leander", "Cheiracanthium_inclusum", "Nitocra_affinis", "Erigone_atra", "Ostrea_edulis")
substitute_taxa<-c("Nabis_rufusculus", "Balaustium_murorum", "Cheiracanthium_punctorium", "Nitocra_hibernica", "Erigone_aletris", "Ostrea_denselamellosa")
for(i in 1:length(missing_taxa)){
	pseudo_data$Consumer.Sp[which(data$Consumer.Sp == missing_taxa[i])]<-substitute_taxa[i]
}

# Check
cbind(sort(unique(data$Consumer.Sp)), sort(unique(pseudo_data$Consumer.Sp)))
# Looks good

# Get the taxonomic info
unique_species<-sort(unique(data$Consumer.Sp)) 
unique_pseudo<-sort(unique(pseudo_data$Consumer.Sp))
ott_ids<-tnrs_match_names(names=unique_pseudo)$ott_id
info<-taxonomy_taxon_info(ott_ids, include_lineage=T)
lineage<-tax_lineage(info)

# Create a summary taxon table 
taxonomic_data<-data.frame(Kingdom=NA, Phylum=NA, Class=NA, Order=NA, Family=NA, Species=unique_species, n_articles=NA, n_experiments=NA)
for(i in 1:nrow(taxonomic_data)){
  kpcof<-lineage[[i]][match(c("kingdom", "phylum", "class", "order", "family"), lineage[[i]]$rank), 2]
  taxonomic_data[i,c(1:5)]<-kpcof
  # How many experiments for each species
  taxonomic_data[i,7]<-length(unique(data$Article.ID[which(data$Consumer.Sp == taxonomic_data$Species[[i]])]))
  taxonomic_data[i,8]<-length(unique(data$Comparison.ID[which(data$Consumer.Sp == taxonomic_data$Species[[i]])]))
}
write.table(taxonomic_data, file="tables/taxonomic_data.csv", row.names=F, col.names=names(taxonomic_data), sep=",")

# We will make a tree for each subset of the data (i.e., each trait type)
traits<-unique(data$Trait_type)
plot_titles<-c("Development Time", "Lifespan", "Reproductive Function", "Body Size")

# List to hold correlation matrices derived from trees
cov_mats<-list()

# Plot the trees
pdf("plots/figure_S3.pdf", height=10, width=10)
par(mfrow=c(2,2), mar=c(0.5, 0.5, 1.5, 0.5))

# Loop for each trait type
for(i in 1:length(traits)){

	# Create a datasubset for the trait (using the psuedo dataset with substituted taxa)
	data_i<-pseudo_data[which(pseudo_data$Trait_type == traits[i]),]
	data_i<-droplevels(data_i)
	
	# Get the taxa from the subset
	sp<-unique(data_i$Consumer.Sp)

	# # Get the tree and save in newick format
	taxa<-tnrs_match_names(names=sp)
	tree<-tol_induced_subtree(ott_ids=ott_id(taxa), file=paste0("data/", traits[i], "_tree.tre"))

	# Read in using ape and convert to ultrametric then to cor matrix
	tree<-read.tree(paste0("data/", traits[i], "_tree.tre"))
	tree<-compute.brlen(tree)
	CovMatrix<-vcv(tree, corr=T)
		
	# Now we have the covariance matrix, though some of the species names no longer match the original data so we need to recode them
	
	# Remove ott IDs
	sp<-rownames(CovMatrix)
	for(j in 1:length(sp)){
		split<-strsplit(sp[j], "_")[[1]]
		split<-split[-length(split)]
		sp[j]<-paste0(split, collapse=" ")
	}
	rm(split)
	rm(j)

	# Need to match the retrieved species names from the taxa object back to the original names. 
	# Note this is not the sister species substitution. The tnrs_match_names function sometimes substitutes species names for synonyms
	new_sp<-taxa$search_string[match(sp, taxa$unique_name)]
	# Convert to upper case to match the original dataset
	new_sp<-str_to_title(new_sp)

	# Now check if we have any substituted names, and if so, substitute back
	for(j in 1:length(missing_taxa)){
		tag<-which(new_sp == substitute_taxa[j])
		if(length(tag) > 0){
			new_sp[tag]<-missing_taxa[j]
		}
	}
	rm(tag)
	rm(j)
		
	# Put the new names in the correlation matrix
	rownames(CovMatrix)<-new_sp
	colnames(CovMatrix)<-new_sp
	
	# Double check we have every species in there
	sp_present<-sort(unique(data$Consumer.Sp[which(data$Trait_type == traits[i])]))
	print(which((sp_present == sort(rownames(CovMatrix))) == F))
	
	# Save the matrix in the list
	cov_mats[[i]]<-CovMatrix
	
	# Plot the tree and save for supp matts
	tree$tip.label<-new_sp
  plot(tree, cex=0.5, main=plot_titles[i])
	
	
	# Clean house
	rm(sp)
	rm(new_sp)
	rm(CovMatrix)
	rm(sp_present)
	rm(taxa)
	rm(tree)
	rm(data_i)
}
dev.off()

# Add names for each phylogeny
names(cov_mats)<-traits
 
# Save the CorMatrix
saveRDS(cov_mats, file="data/PhyloMatrix.Rds")

