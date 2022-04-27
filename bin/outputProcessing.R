#!/usr/bin/env Rscript --vanilla

#------------------------------------------
# Extract command-line arguments

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
out <- args[2]
scov <- as.numeric(args[3])
ident <- as.numeric(args[4])

#------------------------------------------
# Load ECTyper results tsv file
ectyper <- read.table(paste0(out, "/serotype/serotype_", name, "/output.tsv"),
                 header = TRUE, sep = "\t"
                 )

#------------------------------------------
# Load mlst results tsv file

mlst <- read.table(paste0(out, "/mlst/mlst_results_", name, ".tsv"),
                 header = FALSE, sep = "\t",
                 col.names = c("name", "scheme", "ST", "adk", "fumC", 
                               "gyrB", "icd", "mdh", "purA", "recA")
                 )

ST <- mlst[, "ST"]

#------------------------------------------
# Load BLAST results tsv file

blast <- read.table(paste0(out, "/blast/blast_results_", name, ".tsv"),
                    header = FALSE, sep = "\t",
                    col.names = c("Sequence", "Gene", 
                                  "GeneLength", "AlignmentLength", "Mismatches", "Gaps", 
                                  "SequenceStart", "SequenceEnd", "GeneStart", "GeneEnd", 
                                  "Identity", "Evalue", "Bitscore")
                    )

# Calculate subject coverage
blast$"Coverage" <- (blast$AlignmentLength - blast$Gaps) / blast$GeneLength * 100

#------------------------------------------
# Create output file of all APEC ref seqs identified
# (using either user-defined subject coverage and identity or default values)

blastFilterAll <- subset(blast, Coverage >= scov & 
                         Identity >= ident
                        )

write.table(blastFilterAll, file = paste0(out, "/blast_results_", name, ".tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE)

#------------------------------------------
# Identify APEC pathotype markers

blastFilterID <- subset(blast, "Coverage" >= 90 & 
                             "Identity" >= 90 & 
                             Gene %in% c("ompTp|plasmid-encoded_outer_membrane_protease|AY545598.5",
                                         "hlyF|avian_hemolysin_and_antimicrobial_peptide_degradation|DQ381420.1")
                           )

markers <- unique(gsub("\\|.*", "", blastFilterID$Gene))

#------------------------------------------
# Verify sample is E. coli and identify serotype

qc <- ectyper[1, "QC"]
species <- ectyper[1, "Species"]

if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
  serotype <- ectyper[1, "Serotype"]
  Otype <- ectyper[1, "O.type"]
  
  #------------------------------------------
  # Check for presence of APEC plasmid
  
  if ("ompTp" %in% markers & "hlyF" %in% markers) {
    plasmid <- "Present"
  }else {
    plasmid <- "Absent"
  }
  
  #------------------------------------------
  # Identify sequence type and check whether "high risk"

  if (!is.numeric(ST)){
    highRiskST <- "Unknown ST"
  }else if (ST %in% c(131, 23, 428, 355)) {
    highRiskST <- "Yes"
  }else if (Otype == "O78") {
    highRiskST <- "Yes"
  }else if (Otype == "-") {
    highRiskST <- "Unknown"  
  }else {
    highRiskST <- "No"
  }
  
  #------------------------------------------
  # Assign pathotype

  if (!is.numeric(ST)){
    pathotype <- "Unknown ST - pathotype could not be determined."
  }else if (plasmid == "Present" & highRiskST == "Yes"){
    pathotype <- "High Risk APEC"
  }else if (plasmid == "Present" & highRiskST == "No"){
    pathotype <- "APEC"
  }else if (plasmid == "Present" & highRiskST == "Unknown"){
    pathotype <- "APEC - unknown if high risk APEC as O antigen typing failed"
  }else if (plasmid == "Absent" & highRiskST == "Yes"){
    pathotype <- "High Risk non-APEC"
  }else if (plasmid == "Absent" & highRiskST == "No"){
    pathotype <- "non-APEC"
  }else if (plasmid == "Absent" & highRiskST == "Unknown"){
    pathotype <- "non-APEC - unknown if high risk non-APEC as O antigen typing failed"
  }

# If sample is NOT E. coli...
}else {
    serotype <- "NA"
    Otype <- "NA"
    plasmid <- "NA"
    pathotype <- "Not E. coli"
}

#------------------------------------------
# Write results to tsv file

df <- data.frame("Sample" = name, "Species" = species,
                 "Serotype" = serotype, "SerotypeQC" = qc,
                 "ST" = ST, 
                 "APEC_plasmid" = plasmid, "Pathotype" = pathotype)

write.table(df, file = paste0(out, "/pathotype_results_", name, ".tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE)
