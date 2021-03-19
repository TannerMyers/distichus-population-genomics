setwd("~/Dropbox/Distichus_Project/ddRADseq_Phylogeography/stacks/info")

popmap <- read.table("../popmap.tsv", header=FALSE, sep="\t")

# Only need lineage ID
popmap <- popmap[,c(1,3)]

levels(popmap[,2])
#[1] "altavelensis" "aurifer"      "biminiensis"  "brevirostris" "caudalis"    
#[6] "dapsilis"     "distichoides" "distichus"    "dom1"         "dom12"       
#[11] "dom2"         "dom3"         "dom4"         "favillarum"   "ignigularis" 
#[16] "igprop"       "marron"       "ocior"        "porcatus"     "properus"    
#[21] "ravig"        "ravitergum"   "sejunctus"    "suppar"       "vinosus"     
#[26] "websteri"


# Write tables for each lineage ID
write.table(popmap[popmap$V3=="favillarum",], "lineage-popmaps/popmap_fav.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="ravitergum",], "lineage-popmaps/popmap_rav.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="ravig",], "lineage-popmaps/popmap_ravig.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="ignigularis",], "lineage-popmaps/popmap_ig.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="properus",], "lineage-popmaps/popmap_prop.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="igprop",], "lineage-popmaps/popmap_igprop.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="sejunctus",], "lineage-popmaps/popmap_sej.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dom1",], "lineage-popmaps/popmap_dom1.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dom2",], "lineage-popmaps/popmap_dom2.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dom12",], "lineage-popmaps/popmap_dom12.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dom3",], "lineage-popmaps/popmap_dom3.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dom4",], "lineage-popmaps/popmap_dom4.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="aurifer",], "lineage-popmaps/popmap_aurifer.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="suppar",], "lineage-popmaps/popmap_suppar.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="vinosus",], "lineage-popmaps/popmap_vinosus.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="brevirostris",], "lineage-popmaps/popmap_brev.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="caudalis",], "lineage-popmaps/popmap_caudalis.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="websteri",], "lineage-popmaps/popmap_websteri.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="marron",], "lineage-popmaps/popmap_marron.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="altavelensis",], "lineage-popmaps/popmap_altavelensis.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="biminiensis",], "lineage-popmaps/popmap_biminiensis.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="dapsilis",], "lineage-popmaps/popmap_dapsilis.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="distichus",], "lineage-popmaps/popmap_distichus.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="distichoides",], "lineage-popmaps/popmap_distichoides.tsv", sep="\t", row.names=F, col.names=F, quote=F)
write.table(popmap[popmap$V3=="ocior",], "lineage-popmaps/popmap_ocior.tsv", sep="\t", row.names=F, col.names=F, quote=F)
