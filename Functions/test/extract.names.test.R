#Testing

#Data loading
matrix<-read.nexus.data("../../Data/Matrices/GL1988-GH.nex")
tree_list<-read.table("../../Data/Taxon_References/FritzTree_taxa.txt", header=F, stringsAsFactors=F)
WR_list<-read.csv("../../Data/Taxon_References/WilsonReederMSW.csv", header=T, stringsAsFactors=F)
Paleo_list<-read.csv("../../Data/Taxon_References/pdb_taxa.csv", header=T, stringsAsFactors=F)
Ref_list<-list(tree_list, WR_list, Paleo_list)
names(Ref_list)<-c("tree", "WR", "Paleo")


#Find match
expect_is(find.match("Bob_the_builder", Ref_list[[1]]), 'logical') ; message('.', appendLF=FALSE)
expect_false(find.match("Bob_the_builder", Ref_list[[1]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Bob_the_builder", Ref_list[[2]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Bob_the_builder", Ref_list[[3]])) ; message('.', appendLF=FALSE)

expect_true(find.match("Tachyglossus_aculeatus", Ref_list[[1]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Tachyglossus_aculeatus", Ref_list[[2]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Tachyglossus_aculeatus", Ref_list[[3]])) ; message('.', appendLF=FALSE)

expect_false(find.match("Suidae", Ref_list[[1]])) ; message('.', appendLF=FALSE)
expect_true(find.match("Suidae", Ref_list[[2]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Suidae", Ref_list[[3]])) ; message('.', appendLF=FALSE)

expect_false(find.match("Amphiperatherium", Ref_list[[1]])) ; message('.', appendLF=FALSE)
expect_false(find.match("Amphiperatherium", Ref_list[[2]])) ; message('.', appendLF=FALSE)
expect_true(find.match("Amphiperatherium", Ref_list[[3]])) ; message('.', appendLF=FALSE)

#level.column
is_match<-grep("Bob_the_builder", as.character(unlist(Ref_list[[1]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(length(level.column(is_match, Ref_list[[1]])), 0) ; message('.', appendLF=FALSE)
is_match<-grep("Bob_the_builder", as.character(unlist(Ref_list[[2]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(length(level.column(is_match, Ref_list[[2]])), 0) ; message('.', appendLF=FALSE)
is_match<-grep("Bob_the_builder", as.character(unlist(Ref_list[[3]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(length(level.column(is_match, Ref_list[[3]])), 0) ; message('.', appendLF=FALSE)

is_match<-grep("Tachyglossus_aculeatus", as.character(unlist(Ref_list[[1]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(level.column(is_match, Ref_list[[1]]), 1) ; message('.', appendLF=FALSE)

is_match<-grep("Suidae", as.character(unlist(Ref_list[[2]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(level.column(is_match, Ref_list[[2]]), 5) ; message('.', appendLF=FALSE)

is_match<-grep("Amphiperatherium", as.character(unlist(Ref_list[[3]])), ignore.case=TRUE, value=TRUE)[1]
expect_equal(level.column(is_match, Ref_list[[3]]), 1) ; message('.', appendLF=FALSE)

#check.split
expect_equal(check.split("Bob@the@builder", Ref_list[[1]]), c(NA, NA))  ; message('.', appendLF=FALSE)
expect_equal(check.split("Bob@the@builder", Ref_list[[1]]), c(NA, NA))  ; message('.', appendLF=FALSE)
expect_equal(check.split("Bob@the@builder", Ref_list[[2]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Bob@the@builder", Ref_list[[3]]), c(NA, NA)) ; message('.', appendLF=FALSE)

expect_equal(check.split("Bob_the_builder", Ref_list[[1]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Bob_the_builder", Ref_list[[2]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Bob_the_builder", Ref_list[[3]]), c(NA, NA)) ; message('.', appendLF=FALSE)

expect_equal(check.split("Tachyglossus_aculeatus", Ref_list[[1]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Tachyglossus_aculeatus", Ref_list[[2]]), c("Tachyglossus", "aculeatus")) ; message('.', appendLF=FALSE)
expect_equal(check.split("Tachyglossus_aculeatus", Ref_list[[3]]), c("Tachyglossus", "aculeatus")) ; message('.', appendLF=FALSE)

expect_equal(check.split("Amphiperatherium_maximum", Ref_list[[1]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Amphiperatherium_maximum", Ref_list[[2]]), c(NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(check.split("Amphiperatherium_maximum", Ref_list[[3]]), c("Amphiperatherium", "maximum")) ; message('.', appendLF=FALSE)

#match taxon
expect_equal(match.taxon("Bob@the@builder", Ref_list), c(NA, NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Bob_the_builder", Ref_list), c(NA, NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Tachyglossus_aculeatus", Ref_list), c("TRUE","tree","Species")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("T._aculeatus", Ref_list), c(NA, NA, NA)) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Tachyglossus", Ref_list), c("TRUE","WR","Genus")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("tACHYGLOSSUS", Ref_list), c("TRUE","WR","Genus")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("aculeatus", Ref_list), c("TRUE","WR","Species")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Suidae", Ref_list), c("TRUE","WR","Family")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Amphiperatherium_maximum", Ref_list), c("FALSE","Paleo","Species")) ; message('.', appendLF=FALSE)
expect_equal(match.taxon("Amphiperatherium", Ref_list), c("FALSE","Paleo","Genus")) ; message('.', appendLF=FALSE)

#extract name
matching_table<-extract.names(matrix, Ref_list)
expect_is(matching_table, "data.frame") ; message('.', appendLF=FALSE)