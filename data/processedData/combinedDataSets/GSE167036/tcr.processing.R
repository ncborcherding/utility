getwd()
dat <- read.csv("GSE167036_meta_tcr.csv")

CTnt <- paste0(dat$TRA_1_cdr3_nt, "_", dat$TRB_1_cdr3_nt)
CTaa <- paste0(dat$TRA_1_cdr3, "_", dat$TRB_1_cdr3)
CTgene <- paste0(dat$TRA_1_v_gene, ".", dat$TRA_1_j_gene, ".", dat$TRA_1_c_gene, "_",
                 dat$TRB_1_v_gene, ".", dat$TRB_1_d_gene, ".", dat$TRB_1_j_gene, ".", dat$TRB_1_c_gene)
CTstrict <- paste0(dat$TRA_1_v_gene, ".", dat$TRA_1_j_gene, ".", dat$TRA_1_c_gene, "_", dat$TRA_1_cdr3_nt, "_",
                  dat$TRB_1_v_gene, ".", dat$TRB_1_d_gene, ".", dat$TRB_1_j_gene, ".", dat$TRB_1_c_gene, "_", dat$TRB_1_cdr3_nt)

processed.data <- data.frame(CTgene, CTnt, CTaa, CTstrict, sample = dat$patient_id, ID = substring(dat$orig.ident,1,2))
processed.data[processed.data == "nan_nan"] <- NA
processed.data[which(is.na(processed.data$CTnt)), "CTgene"] <- NA
processed.data[which(is.na(processed.data$CTnt)), "CTstrict"] <- NA
processed.data$barcode <- dat$cell_id


processed.data <- processed.data %>%
    group_by(sample, CTnt) %>%
    mutate(Frequency = n()) %>%
    as.data.frame()



processed.data <- na.omit(processed.data)
rownames(processed.data) <- processed.data$barcode

nsize <- table(processed.data$sample)
uniq.patients <- unique(processed.data$sample)

for(i in seq_along(uniq.patients)) {
    processed.data[processed.data$sample == uniq.patients[i],"Frequency"] <- processed.data[processed.data$sample == uniq.patients[i],"Frequency"]/nsize[i]
}

cloneTypes=c(None = 0, Rare = 1e-4, Small = 0.001, 
             Medium = 0.01, Large = 0.1, Hyperexpanded = 1)


processed.data$cloneType <- NA
for (x in seq_along(cloneTypes)) { names(cloneTypes)[x] <- 
    paste0(names(cloneTypes[x]), ' (', cloneTypes[x-1], 
           ' < X <= ', cloneTypes[x], ')') }



for (i in 2:length(cloneTypes)) { processed.data$cloneType <- 
    ifelse(processed.data$Frequency > cloneTypes[i-1] & processed.data$Frequency 
           <= cloneTypes[i], names(cloneTypes[i]), processed.data$cloneType) }

saveRDS(processed.data, "ProcessedTCR_forutility.rds")

