# fisher test
df.fisher <- data.frame()
feqtl.raqtl <- '/local/zy/PEI/raQTL/liver/filterd_eqtl_raqtl.txt'
df.feqtl.raqtl <- read.delim(feqtl.raqtl, sep = '\t', stringsAsFactors = F)
df.input <- df.feqtl.raqtl[, c("snp_eqtl", "snp_raqtl")]
df.input[df.input$snp_eqtl != '.', "snp_eqtl"] = 'eQTL'
df.input[df.input$snp_eqtl == '.', "snp_eqtl"] = 'non-eQTL'
df.input[df.input$snp_raqtl != '.',"snp_raqtl"] = 'raQTL'
df.input[df.input$snp_raqtl == '.',"snp_raqtl"] = 'non-raQTL'
table.input <- xtabs(~ snp_eqtl+snp_raqtl, data = df.input)
