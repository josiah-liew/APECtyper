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

qc <- ectyper[1, "QC"]
species <- ectyper[1, "Species"]

#------------------------------------------
# Load mlst results tsv file

mlst <- read.table(paste0(out, "/mlst/mlst_results_", name, ".tsv"),
                 header = FALSE, sep = "\t"
                 )

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

if (grepl("Escherichia coli", species, ignore.case = TRUE)) {
  serotype <- ectyper[1, "Serotype"]
  Otype <- ectyper[1, "O.type"]
  ST <- mlst[, 3]
  
  #------------------------------------------
  # Check for presence of APEC plasmid
  
  if ("ompTp" %in% markers & "hlyF" %in% markers) {
    plasmid <- "Present"
  }else {
    plasmid <- "Absent"
  }
  
  #------------------------------------------
  # Check whether "high risk"
  
  if (Otype == "O78" || ST %in% c(131, 23, 428, 355)) {
    highRisk <- "Yes"
  }else if (Otype == "-" && !is.numeric(ST)) {
    highRisk <- "UnknownSTO"
  }else if (Otype == "-" && is.numeric(ST)) {
    highRisk <- "UnknownO"
  }else if (!is.numeric(ST)) {
    highRisk <- "UnknownST"
  }else if (is.numeric(ST)) {
    highRisk <- "No"
  }
  
  #------------------------------------------
  # Assign pathotype

  if (plasmid == "Present") {
    if (highRisk == "Yes") {
      pathotype <- "High Risk APEC"
    }else if (highRisk == "No") {
      pathotype <- "APEC"
    }else if (highRisk == "UnknownSTO") {
      pathotype <- "APEC - unknown whether isolate is high risk (unknown ST and O-type)"
    }else if (highRisk == "UnknownO") {
      pathotype <- "APEC - unknown whether isolate is high risk (unknown O-type)"       
    }else if (highRisk == "UnknownST") {
      pathotype <- "APEC - unknown whether isolate is high risk (unknown ST)"
    }
  }else if (plasmid == "Absent") {
    if (highRisk == "Yes") {
      pathotype <- "High Risk non-APEC"
    }else if (highRisk == "No") {
      pathotype <- "non-APEC"
    }else if (highRisk == "UnknownSTO") {
      pathotype <- "non-APEC - unknown whether isolate is high risk (unknown ST and O-type)"
    }else if (highRisk == "UnknownO") {
      pathotype <- "non-APEC - unknown whether isolate is high risk (unknown O-type)"       
    }else if (highRisk == "UnknownST") {
      pathotype <- "non-APEC - unknown whether isolate is high risk (unknown ST)"
    }
  }

# If sample is NOT E. coli...
}else {
  serotype <- "NA"
  Otype <- "NA"
  ST <- "NA"
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
