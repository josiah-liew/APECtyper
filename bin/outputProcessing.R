#!/usr/bin/env Rscript --vanilla

#------------------------------------------
# Extract command-line arguments
args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
out <- args[2]
scov <- as.numeric(args[3])
ident <- as.numeric(args[4])

name
out
scov
ident

# name="PP394_S59"
# out=
# scov=90
# ident=90

#------------------------------------------

# mlst <- read.csv("/Users/elizabethmiller/Desktop/Projects/APECtyper/out_multiSeq/mlst/mlst_results_PP394_S59.csv",
#                  header = FALSE, 
#                  col.names = c("name", "scheme", "ST", "adk", "fumC", 
#                                "gyrB", "icd", "mdh", "purA", "recA")
# )


# Load mlst results csv file
mlst <- read.csv(paste0(out, "/mlst/mlst_results_", name, ".csv"),
                 header = FALSE, 
                 col.names = c("name", "scheme", "ST", "adk", "fumC", 
                               "gyrB", "icd", "mdh", "purA", "recA")
                 )

# Check that mlst file has exactly one row
if (nrow(mlst) != 1) {
  print("Error in mlst output.")
}

mlst

ST <- mlst[, "ST"]

ST

#------------------------------------------
# Load BLAST results tsv file

# blast <- read.table("/Users/elizabethmiller/Desktop/Projects/APECtyper/out_multiSeq/blast/blast_results_PP394_S59.tsv",
#                     header = FALSE, sep = "\t", 
#                     col.names = c("Sequence", "Gene", 
#                                   "GeneLength", "AlignmentLength", "Mismatches", "Gaps", 
#                                   "SequenceStart", "SequenceEnd", "GeneStart", "GeneEnd", 
#                                   "Identity", "Evalue", "Bitscore")
#                     )

blast <- read.table(paste0(out, "/blast/blast_results_", name, ".tsv"),
                    header = FALSE, sep = "\t",
                    col.names = c("Sequence", "Gene", 
                                  "GeneLength", "AlignmentLength", "Mismatches", "Gaps", 
                                  "SequenceStart", "SequenceEnd", "GeneStart", "GeneEnd", 
                                  "Identity", "Evalue", "Bitscore")
                    )

# Calculate subject coverage
blast$"Coverage" <- (blast$AlignmentLength - blast$Gaps) / blast$GeneLength * 100
blast

#------------------------------------------
# Create output file of all APEC ref seqs identified
# (using either user-defined subject coverage and identity or default values)

blastFilterAll <- subset(blast, Coverage >= scov & 
                      Identity >= ident & 
                      Gene != "O78|O-antigen_gene_cluster_partial|FJ940775.1"
                    )

write.csv(blastFilterAll, file = paste0(out, "/blast_results_", name, ".csv"),
          row.names = FALSE, quote = FALSE)

#------------------------------------------
# Identify APEC markers and define APEC pathotype

blastFilterID <- subset(blast, "Coverage" >= 90 & 
                             "Identity" >= 90 & 
                             Gene %in% c("O78|O-antigen_gene_cluster_partial|FJ940775.1",
                                         "ompTp|plasmid-encoded_outer_membrane_protease|AY545598.5",
                                         "hlyF|avian_hemolysin_and_antimicrobial_peptide_degradation|DQ381420.1")
                           )

markers <- unique(gsub("\\|.*", "", blastFilterID$Gene))
markers

if ("ompTp" %in% markers & "hlyF" %in% markers) {
  plasmid <- "Present"
  }else {
    plasmid <- "Absent"
    }

paste0("APEC plasmid: ", plasmid)

if ("O78" %in% markers) {
  O78 <- "O78"
  }else {
    O78 <- "Not O78"
    }

paste0("Serogroup: ", O78)

if (ST %in% c(131, 23, 428, 355) || O78 == "Found" || ST == 117 & O78 == "Found") {
  highRiskST <- "Present"
  }else {
    highRiskST <- "Absent"
    }

paste0("High Risk ST: ", highRiskST)


if (plasmid == "Present" & highRiskST == "Present"){
  pathotype <- "High Risk APEC"
  }else if (plasmid == "Present" & highRiskST == "Absent"){
    pathotype <- "APEC"
    }else if (plasmid == "Absent" & highRiskST == "Present"){
      pathotype <- "High Risk non-APEC"
      }else if (plasmid == "Absent" & highRiskST == "Absent"){
        pathotype <- "non-APEC"
        }

paste0("Pathotype: ", pathotype)

ST
Serogroup
hlyF
ompT
Pathotype

df <- data.frame(Sample = name, "ST" = ST, "Serogroup" = O78, "APEC plasmid" = plasmid, "Pathotype" = pathotype)

write.csv(df, file = paste0(out, "/pathotype_results_", name, ".csv"),
          row.names = FALSE, quote = FALSE)
