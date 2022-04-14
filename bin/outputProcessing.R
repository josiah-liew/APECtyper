#!/usr/bin/env Rscript --vanilla

#------------------------------------------
# Extract command-line arguments

args <- commandArgs(trailingOnly = TRUE)
name <- args[1]
out <- args[2]
scov <- as.numeric(args[3])
ident <- as.numeric(args[4])

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
                      Identity >= ident & 
                      Gene != "O78|O-antigen_gene_cluster_partial|FJ940775.1"
                    )

write.table(blastFilterAll, file = paste0(out, "/blast_results_", name, ".tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE)

#------------------------------------------
# Identify APEC markers and define APEC pathotype

blastFilterID <- subset(blast, "Coverage" >= 90 & 
                             "Identity" >= 90 & 
                             Gene %in% c("O78|O-antigen_gene_cluster_partial|FJ940775.1",
                                         "ompTp|plasmid-encoded_outer_membrane_protease|AY545598.5",
                                         "hlyF|avian_hemolysin_and_antimicrobial_peptide_degradation|DQ381420.1")
                           )

markers <- unique(gsub("\\|.*", "", blastFilterID$Gene))

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

if (!is.numeric(ST)){
  highRiskST <- "Error"
  }else if (ST %in% c(131, 23, 428, 355) || O78 == "Found") {
    highRiskST <- "Present"
    }else {
      highRiskST <- "Absent"
      }

paste0("High Risk ST: ", highRiskST)

if (!is.numeric(ST)){
  pathotype <- "ST could not be identified - check if isolate is really E. coli"
  }else if (plasmid == "Present" & highRiskST == "Present"){
    pathotype <- "High Risk APEC"
    }else if (plasmid == "Present" & highRiskST == "Absent"){
      pathotype <- "APEC"
      }else if (plasmid == "Absent" & highRiskST == "Present"){
        pathotype <- "High Risk non-APEC"
        }else if (plasmid == "Absent" & highRiskST == "Absent"){
          pathotype <- "non-APEC"
          }

paste0("Pathotype: ", pathotype)

df <- data.frame(Sample = name, "ST" = ST, "Serogroup" = O78, "APEC plasmid" = plasmid, "Pathotype" = pathotype)

write.table(df, file = paste0(out, "/pathotype_results_", name, ".tsv"),
          sep = "\t", row.names = FALSE, quote = FALSE)
