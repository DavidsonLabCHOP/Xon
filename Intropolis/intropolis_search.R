# This script is released as a part of Monteys et al., Nature (2021)
# It is used to query the presence of a splice junction of interest in the Intropolis database of exon-exon junctions.
# Intropolis contains a list of exon-exon junctions found across 21,504 RNA-Seq Datasets
# It can be accessed here https://github.com/nellore/intropolis
# The purpose of this script is to identify lines in the intropolis database corresponding to a splicing event of interest 
# To do this you need to have access to a copy of the intropolis database containing only the first 4 columns (rowNumber, chromosome, startPosition, endPosition).
# Once the following script is used to identify the line of the intropolis database represinting a splice event of interest,
# we then use a second python-based script "ExtractIntropolisLine.py" to extract (by rowNumber) the number and identity of datasets that contained that splice event.
library(readr)

# To begin download the large intropolis database "intropolis.v1.hg19.tsv.gz" and uncompress. 
# Here we name this uncompressed database "smallIntropDB.tsv" and store it in our downloads folder on our mac machine.
smallIntropDB <- data.frame(read_table2("~/Downloads/smallIntropDB.tsv"), stringsAsFactors = FALSE)
smallIntropDB <- data.frame(smallIntropDB[,1], smallIntropDB[,2], smallIntropDB[,3], smallIntropDB[,4])
colnames(smallIntropDB) <- c("chr","start","end","NA")
head(smallIntropDB)

# Sf3b3 Genomic Minigene Region
# chr16:70,526,657-70,529,199
# chr5:55243418-55247291

chr <- c("chr16")
start <- c(70526657-100:70526657+100)
end <- c(70526657-100:70529199+100)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetChrStart <- subsetQuery[subsetQuery$start %in% start,]
subsetChrEnd <- subsetQuery[subsetQuery$end %in% end,]

# Sf3b3 Pseudo exon position
# chr16:70,527,376-70,527,429

chr <- c("chr16")
start <- c(70526657:70527500)
end <- c(70528500:70529199)

subsetQuery <- data.frame(smallIntropDB[smallIntropDB$chr == chr,], stringsAsFactors = FALSE)
subsetStartPos <- subsetQuery[subsetQuery$start %in% c(70526657:70527500),]
subsetEndPos <- subsetQuery[subsetQuery$end %in% end,]


# Test plot generation
chr <- c("chr16")
start <- c(70526689:70526774)
end <- c(70528831:70528911)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetStartPos <- subsetQuery[subsetQuery$start %in% start,]
subsetEndPos <- subsetStartPos[subsetStartPos$end %in% end,]


# ACTB example
chr <- c("chr7")
GenePos <- c(5527148:5530601)
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% GenePos,]

# ACTB Exon1 example
chr <- c("chr7")
start <- c(5527870:5527911)
end <- c(5527995:5528012)
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetQuery[subsetQuery$start %in% end,]

# SF3B3 example
chr <- c("chr16")
GenePos <- c(70523816:70577668)
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% GenePos,]

##############
# SF3B3
##############
# SF3B3 5' intron of novel exon example
chr <- c("chr16")
start <- c(70560625:70560635)
end <- c(70561274:70561284)
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]

# SF3B3 3' intron of novel exon example
chr <- c("chr16")
start <- c(c(70561332 - 5):c(70561332 + 5))
end <- c(70561279:70561332)
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]

################
# BENC1
################
# BENC1 5' intron of novel exon example
chr <- "chr12"
pos1 <- 42481737 # The 5' most splice donor
pos2 <- 42488953 # The 5' novel exon edge
pos3 <- 42489016 # The 3' novel exon edge
pos4 <- 42481588 # The 3' most splice donor

start <- c(c(pos1 - 5):c(pos1 + 5))
end <- c(c(pos2 - 5):c(pos2 + 5))
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]
exon1_result

# BENC1 3' intron of novel exon example
start <- c(c(pos3 - 5):c(pos3 + 5))
end <- c(c(pos4 - 5):c(pos4 + 5))
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]
exon1_result

################
# Actb
################
# Actb 5' intron of novel exon example
pos1 <- 5567522 # The 5' most splice donor
pos2 <- 5567634 # The 5' novel exon edge
pos3 <- 5567816 # The 3' novel exon edge
pos4 <- 5567911 # The 3' most splice donor

chr <- c("chr7")
start <- c(c(pos1 - 5):c(pos1 + 5))
end <- c(c(pos2 - 5):c(pos2 + 5))
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]
exon1_result

# BENC1 3' intron of novel exon example
chr <- c("chr17")
start <- c(c(pos3 - 5):c(pos3 + 5))
end <- c(c(pos4 - 5):c(pos4 + 5))
#end <- c(5528000:5528009)

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
exon1_result <- subsetGenePos[subsetGenePos$end %in% end,]
exon1_result


#
##
###
####
####################
# Search Intropolis function
####################

search_intropolis_fun <- function(chr, pos1, pos2, pos3, pos4) {

start <- c(c(pos1 - 5):c(pos1 + 5))
end <- c(c(pos2 - 5):c(pos2 + 5))

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
intron1_result <- subsetGenePos[subsetGenePos$end %in% end,]
intron1_result

start <- c(c(pos3 - 5):c(pos3 + 5))
end <- c(c(pos4 - 5):c(pos4 + 5))

subsetQuery <- smallIntropDB[smallIntropDB$chr == chr,]
subsetGenePos <- subsetQuery[subsetQuery$start %in% start,]
intron2_result <- subsetGenePos[subsetGenePos$end %in% end,]
intron2_result
result <- list(intron1_result, intron2_result)
return(result)
}

Actb_example1 <- search_intropolis_fun("chr7", 5568350, 5568649, 5568690, 5568791)
Actb_example2 <- search_intropolis_fun("chr7", 5568350, 5568791, 5569030, 5569165)

SF3B3_std_result <- search_intropolis_fun("chr16", 70560629, 70562775, 70560629, 70562775)
SF3B3_result <- search_intropolis_fun("chr16", 70560629, 70561279, 70561332, 70562775)
BENC1_result <- search_intropolis_fun("chr17", 40962946, 40963310, 40963348, 40963672)
GXYLT1_result <- search_intropolis_fun("chr12", 42481749, 42488953, 42489016, 42491243)
SKP1_1_result <- search_intropolis_fun("chr5", 133509716, 133510975, 133511076, 133512545)
SKP1_2_result <- search_intropolis_fun("chr5", 133509716, 133510975, 133511114, 133512545)
C12orf4_result <- search_intropolis_fun("chr12", 4645385, 4646546, 4646680, 4647574)

Myo7a <- search_intropolis_fun("chr11", 76877198, 76883794, 76877198, 76883794)
Myo7a <- search_intropolis_fun("chr11", 76877198, 76877985, 76878038, 76883794)
                      

#
##
###
####################
# Run for spreadsheet
####################
SF3B3_Canonical <- search_intropolis_fun("chr16", 70560629, 70562775, 70560629, 70562775)
BENC1_Canonical <- search_intropolis_fun("chr17", 40962946, 40963672, 40962946, 40963672)
GXYLT1_Canonical <- search_intropolis_fun("chr12", 42481749, 42491243, 42481749, 42491243)
SKP1_Canonical <- search_intropolis_fun("chr5", 133509713, 133512544, 133509713, 133512544)
C12ORF4_Canonical <- search_intropolis_fun("chr12", 4645385, 4647575, 4645385, 4647575)
SSBP1_Canonical <- search_intropolis_fun("chr7", 141438990, 141441968, 141438990, 141441968)
RARS_Canonical <- search_intropolis_fun("chr5", 167945067, 167946085, 167945067, 167946085)
PXDC2P_Canonical <-  search_intropolis_fun("chr16", 70064969, 70065802, 70064969, 70065802)
STRADB_Canonical <- search_intropolis_fun("chr2", 202334775, 202337677, 202334775, 202337677)
WNK1_Canonical <- search_intropolis_fun("chr12", 1003801, 1005236, 1003801, 1005236)
WDR27_Canonical <-  search_intropolis_fun("chr6", 170060862, 170062399, 170060862, 170062399)
Cip2a_Canonical <- search_intropolis_fun("chr3", 108284301, 108285343, 108284301, 108285343)
IFT57_Canonical <- search_intropolis_fun("chr3", 107910490, 107925474, 107910490, 107925474)
HTT_Canonical <- search_intropolis_fun("chr4", 3214436, 3215684, 3214436, 3215684)
SAK2_Canonical <- search_intropolis_fun("chr17", 57189706, 57196679, 57189706, 57196679)
EVC_Canonical <- search_intropolis_fun("chr4", 5735162, 5743442, 5735162, 5743442)
DYRK1A_Canonical <- search_intropolis_fun("chr21", 38792686, 38844985, 38792686, 38844985)
GNAQ_Canonical <- search_intropolis_fun("chr9", 80430686, 80537076, 80430686, 80537076)
ZMYM6_Canonical <- search_intropolis_fun("chr1", 35485203, 35485983, 35485203, 35485983)
CYB5B_Canonical  <- search_intropolis_fun("chr16", 69482047, 69492995, 69482047, 69492995)
MMS22L_Canonical <- search_intropolis_fun("chr6", 97634566, 97676769, 97634566, 97676769)
MEMO1_Canonical <- search_intropolis_fun("chr2", 32108531, 32117060, 32108531, 32117060)
PNISR_Canonical <- search_intropolis_fun("chr6", 99864304, 99873090, 99864304, 99873090)

CACNA2D1_Canonical <- search_intropolis_fun("chr7", 81695840, 81714084, 81695840, 81714084)
SSBP1_Canonical <- search_intropolis_fun("chr7", 141438990, 141441968, 141438990, 141441968)
DDX42_Canonical <- search_intropolis_fun("chr17", 61882535, 61883894, 61882535, 61883894)
ASAP1_Canonical <- search_intropolis_fun("chr8", 131172210, 131179781, 131172210, 131179781)
DUXAP10_Canonical <- search_intropolis_fun("chr14", 19884029, 19894699, 19884029, 19894699)
AVL9_Canonical <- search_intropolis_fun("chr7", 32599074, 32609631, 32599074, 32609631)
DYRK1A_Canonical <- search_intropolis_fun("chr21", 38792686, 38844985, 38792686, 38844985)
FAM3A_Canonical <- search_intropolis_fun("chrX", 153740735, 153741146, 153740735, 153741146)
FHOD3_Canonical <- search_intropolis_fun("chr18", 34320801, 34322699, 34320801, 34322699)
TBCA_Canonical <- search_intropolis_fun("chr5", 77004172, 77072028, 77004172, 77072028)
MZT1_Canonical <- search_intropolis_fun("chr13", 73293235, 73301661, 73293235, 73301661)
LINC_Canonical <- search_intropolis_fun("chr14", 19680685, 19683691, 19680685, 19683691)
SF3B3_2_Canonical <- search_intropolis_fun("chr16", 70575737, 70578340, 70575737, 70578340)
SAFB_Canonical <-  search_intropolis_fun("chr19", 5654211, 5654378, 5654211, 5654378)
GCFC2_Canonical <- search_intropolis_fun("chr2", 75929549, 75933648, 75929549, 75933648)
MRPL45_Canonical <- search_intropolis_fun("chr17", 36462598, 36474585, 36462598, 36474585)
SPIDR_Canonical <- search_intropolis_fun("chr8", 48173583, 48192449,48173583, 48192449)
DUXAP8_Canonical <- search_intropolis_fun("chr14", 19680685, 19683691, 19680685, 19683691)
PDXDC1_Canonical <- search_intropolis_fun("chr16", 15102704, 15103537, 15102704, 15103537)
MAN1A2_Canonical <- search_intropolis_fun("chr1", 117984947, 118003110, 117984947, 118003110)
RAF1_Canonical <- search_intropolis_fun("chr3", 12641914, 12645634, 12641914, 12645634)
ERGIC3_Canonical <- search_intropolis_fun("chr20", 34136617, 34142142, 34136617, 34142142)
IL6ST_Canonnical <- search_intropolis_fun("chr5", 55243418, 55247291, 55243418, 55247291)

chr5:55243418-55247291

#
##
###
####################
# Run Novel Junctions
####################
SF3B3_Novel <- search_intropolis_fun("chr16", 70560629, 70561279, 70561332, 70562775)
BENC1_Novel <- search_intropolis_fun("chr17", 40962946, 40963310, 40963348, 40963672)
GXYLT1_Novel <- search_intropolis_fun("chr12", 42481749, 42488953, 42489016, 42491243)
SKP1_Novel <- search_intropolis_fun("chr5", 133509713, 133510975, 133511076, 133512544)
SKP1_Novel2 <- search_intropolis_fun("chr5", 133509713, 133510975, 133511114, 133512544)
C12ORF4_Novel <- search_intropolis_fun("chr12", 4645385, 4646546, 4646680, 4647575)
SSBP1_Novel <- search_intropolis_fun("chr7", 141438990, 141441110, 141441259, 141441968)
RARS_Novel <- search_intropolis_fun("chr5", 167945067, 167945374, 167945528, 167946085)
RARS_Novel2 <- search_intropolis_fun("chr5", 167945067, 167945474, 167945528, 167946085)

PXDC2P_Novel <-  search_intropolis_fun("chr16", 70064969, 70065089, 70065151, 70065802)
STRADB_Novel <- search_intropolis_fun("chr2", 202334775, 202335630, 202335834, 202337677)
WNK1_Novel <- search_intropolis_fun("chr12", 1003801, 1004327, 1004362, 1005236)
WDR27_Novel <-  search_intropolis_fun("chr6", 170060862, 170061799, 170061846, 170062399)
Cip2a_Novel <- search_intropolis_fun("chr3", 108284301, 108284745, 108284778, 108285343)
IFT57_Novel <- search_intropolis_fun("chr3", 107910490, 107911323, 107911373, 107925474)
HTT_Novel <- search_intropolis_fun("chr4", 3214436, 3215349, 3215463, 3215684)
SAK2_Novel <- search_intropolis_fun("chr17", 57189706, 57196756, 57196856, 57196679)
EVC_Novel <- search_intropolis_fun("chr4", 5735162,5743061,5743168, 5743442)
DYRK1A_Novel <- search_intropolis_fun("chr21", 38792686, 38794883, 38794954, 38844985)
GNAQ_Novel <- search_intropolis_fun("chr9", 80430686, 80535564, 80535619, 80537076)
ZMYM6_Novel <- search_intropolis_fun("chr1", 35485203, 35485862, 35485880, 35485983)
CYB5B_Novel  <- search_intropolis_fun("chr16", 69482047, 69482508, 69482656, 69492995)
MMS22L_Novel <- search_intropolis_fun("chr6", 97634566, 97649238, 97649341, 97676769)
MEMO1_Novel <- search_intropolis_fun("chr2", 32108531, 32112104, 32112156, 32117060)
PNISR_Novel <- search_intropolis_fun("chr6", 99864304, 99868399, 99868460, 99873090)

CACNA2D1_Novel <- search_intropolis_fun("chr7", 81695840, 81705332, 81705438, 81714084)
SSBP1_Novel <- search_intropolis_fun("chr7", 141438990, 141441110, 141441259, 141441968)
DDX42_Novel <- search_intropolis_fun("chr17", 61882535, 61883894, 61882535, 61883894)
ASAP1_Novel <- search_intropolis_fun("chr8", 131172210, 131173031, 131173039, 131179781)
DUXAP10_Novel <- search_intropolis_fun("chr14", 19884029, 19893035, 19893150, 19894699)
AVL9_Novel <- search_intropolis_fun("chr7", 32599074, 32602170, 32602525, 32609631)
DYRK1A_Novel <- search_intropolis_fun("chr21", 38792686, 38794884, 38794954, 38844985)
FAM3A_Novel <- search_intropolis_fun("chrX", 153740735, 153740892, 153741030, 153741146)
FHOD3_Novel <- search_intropolis_fun("chr18", 34320801, 34322340, 34322431, 34322699)
TBCA_Novel <- search_intropolis_fun("chr5", 77004172, 77070041, 77070041, 77072028)
MZT1_Novel <- search_intropolis_fun("chr13", 73293235, 73299780, 73299916, 73301661)
LINC_Novel <- search_intropolis_fun("chr14", 19680685, 19682237, 19682352, 19683691)
SF3B3_2_Novel <- search_intropolis_fun("chr16", 70575737, 70578072, 70578152, 70578340)
SAFB_Novel <-  search_intropolis_fun("chr19", 5654211, 5654151, 5654379, 5654378)
GCFC2_Novel <- search_intropolis_fun("chr2", 75929549, 75929817, 75929933, 75933648)
MRPL45_Novel <- search_intropolis_fun("chr17", 36462598, 36468550, 36468624, 36474585)
SPIDR_Novel <- search_intropolis_fun("chr8", 48173583, 48185929, 48186042, 48192449)
DUXAP8_Novel <- search_intropolis_fun("chr14", 19680685, 19682237, 19682352, 19683691)
PDXDC1_Novel <- search_intropolis_fun("chr16", 15102704, 15103356, 15103418, 15103537)
MAN1A2_Novel <- search_intropolis_fun("chr1", 117984947, 117998707, 117998828, 118003110)
RAF1_Novel <- search_intropolis_fun("chr3", 12641914, 12644977, 12645036, 12645634)
ERGIC3_Novel <- search_intropolis_fun("chr20", 34136617, 34136917, 34136961, 34142142)
IL6ST_Novel <- search_intropolis_fun("chr5", 55243418, 55243929, 55244811, 55247291)
