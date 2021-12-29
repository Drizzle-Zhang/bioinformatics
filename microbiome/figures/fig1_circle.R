library(RColorBrewer)
abundance <- 0.0001
type.cutoff <- 'diff'
topN <- 50
series.time <- c(-1, 1, 5, 9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
pathout <- '/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan'

list.vec.dose <- list(c(0, 1, 2, 3), c(0, 1), c(0, 2), c(0, 3))
vec.gender <- c('male', 'female')
for (vec.dose in list.vec.dose) {
    # dose
    str_dose <- paste0(as.character(vec.dose), collapse = '')
    path.dose <- paste(pathout, str_dose, sep = '/')
    
    file.tree <- 
        paste0(path.dose, '/tree_backbone_', topN, '.txt')
    df.tree <- read.table(file.tree, sep = '.')
    colnames(df.tree) <- c('phylum', 'class', 'order', 'family', 'genus', 'species', 'OTU_ID')
    OTUs.diff <- df.tree$OTU_ID
    
    for (gender in vec.gender) {
        path.diff <- 
            paste0(path.dose, '/', gender, '_', type.cutoff, '_', abundance)
        # set color for each phylum
        table.phylum <- sort(table(df.tree[, 'phylum']), decreasing = T)
        phylums <- names(table.phylum)
        df.color <- data.frame(stringsAsFactors = F)
        # set.color <- brewer.pal(8, "Set2")
        # set.color <- c("#A6CEE3", "#B2DF8A", "#FB9A99", "#FDBF6F", 
        #                "#CAB2D6", "#FFFF99")
        set.color <- c("#6A3D9A", "#FF7F00", "#33A02C", "#E31A1C", 
                       "#1F78B4", "#FFFF99")
        for (i in 1:length(phylums)) {
            if (i > 6 | table.phylum[phylums[i]] < 3) {
                sub.color <- '#D3D3D3'
            } else {
                sub.color <- set.color[i]
            }
            df.color <- rbind(df.color, 
                              data.frame(phylum = phylums[i], color = sub.color,
                                         stringsAsFactors = F))
        }
        df.color <- merge(unique(df.tree[, c('phylum', 'genus')]),
                          df.color, by = 'phylum')
        
        # generate phylum annotation file
        phylum_color <- data.frame(stringsAsFactors = F)
        for (sub.phylum in phylums) {
            sub.anno <- data.frame(stringsAsFactors = F)
            
            sub.anno[1, 1] <- sub.phylum
            sub.anno[1, 2] <- 'annotation_background_color'
            sub.anno[1, 3] <- unique(df.color[df.color$phylum == sub.phylum, 'color'])
            
            sub.anno[2, 1] <- paste0(sub.phylum, '*')
            sub.anno[2, 2] <- 'clade_marker_color'
            sub.anno[2, 3] <- unique(df.color[df.color$phylum == sub.phylum, 'color'])
            
            phylum_color <- rbind(phylum_color, sub.anno)
            
        }
        
        # generate genus annotation file
        families <- unique(df.tree[, 'genus'])
        genus_color <- data.frame(stringsAsFactors = F)
        for (sub.genus in families) {
            # name.genus <- substr(sub.genus, 4, nchar(sub.genus))
            sub.anno <- data.frame(stringsAsFactors = F)
            sub.anno[1, 1] <- sub.genus
            sub.anno[1, 2] <- 'annotation'
            sub.anno[1, 3] <- '*'
            # rotate text 90 degree
            sub.anno[2, 1] <- sub.genus
            sub.anno[2, 2] <- 'annotation_rotation'
            sub.anno[2, 3] <- '90'
            # set background color
            sub.anno[3, 1] <- sub.genus
            sub.anno[3, 2] <- 'annotation_background_color'
            sub.anno[3, 3] <- df.color[df.color$genus == sub.genus, 'color']
            # set font size
            sub.anno[4, 1] <- sub.genus
            sub.anno[4, 2] <- 'annotation_font_size'
            sub.anno[4, 3] <- '15'
            
            genus_color <- rbind(genus_color, sub.anno)
            
        }
        
        label_color <- rbind(phylum_color, genus_color)
        write.table(label_color, paste0(path.dose, '/tree_label_color_', topN, '.txt'), 
                    sep = '\t', quote = F, row.names = F, col.names = F, na = '')
        
        # generate ring annotation file
        limit <- 3
        # color.pos <- '#F8766D'
        # color.neg <- '#00BFC4'
        # color.pos <- '#EE4000'
        # color.neg <- '#20B2AA'
        color.pos <- '#8B0000'
        color.neg <- '#00008B'
        df.ring <- data.frame(stringsAsFactors = F)
        df.ring.global <- data.frame(stringsAsFactors = F)
        for (i in 1:length(series.time)) {
            sub.time <- series.time[i]
            
            time.anno <- data.frame(stringsAsFactors = F)
            time.anno[1, 1] <- 'ring_internal_separator_thickness'
            time.anno[1, 2] <- i
            time.anno[1, 3] <- 0.2
            time.anno[2, 1] <- 'ring_width'
            time.anno[2, 2] <- i
            time.anno[2, 3] <- 1
            time.anno[3, 1] <- 'ring_height'
            time.anno[3, 2] <- i
            time.anno[3, 3] <- 1
            time.anno[4, 1] <- 'ring_label'
            time.anno[4, 2] <- i
            time.anno[4, 3] <- sub.time
            time.anno[5, 1] <- 'ring_label_font_size'
            time.anno[5, 2] <- i
            time.anno[5, 3] <- 17
            time.anno[6, 1] <- 'ring_separator_color'
            time.anno[6, 2] <- i
            time.anno[6, 3] <- 'gray'
            df.ring.global <- rbind(df.ring.global, time.anno)
            
            file.res <- paste0(path.diff, "/OUT_edgeR_", sub.time, 
                               paste0(as.character(vec.dose), collapse = ''), 
                               ".txt")
            res.edgeR <- read.delim(file.res, row.names = 1)
            res.edgeR <- res.edgeR[res.edgeR$Row.names %in% OTUs.diff,]
            res.edgeR$logPval <- log10(res.edgeR$PValue) * 
                (res.edgeR$logFC / abs(res.edgeR$logFC))
            res.edgeR$logPval[is.na(res.edgeR$logPval)] <- 0
            res.edgeR$logPval[abs(res.edgeR$logPval) < 0.5] <- 0
            for (j in 1:dim(res.edgeR)[1]) {
                sub.anno <- data.frame(stringsAsFactors = F)
                OTU_id <- res.edgeR[j, 'Row.names']
                sub.logPval <- res.edgeR[j, "logPval"]
                sub.anno[1, 1] <- OTU_id
                sub.anno[1, 2] <- 'ring_color'
                sub.anno[1, 3] <- i
                if (sub.logPval > 0) {
                    sub.anno[1, 4] <- color.pos
                } else {
                    sub.anno[1, 4] <- color.neg
                }
                sub.anno[2, 1] <- OTU_id
                sub.anno[2, 2] <- 'ring_alpha'
                sub.anno[2, 3] <- i
                sub.anno[2, 4] <- min(1, abs(sub.logPval/limit))
                df.ring <- rbind(df.ring, sub.anno)
            }
        }
        
        write.table(df.ring.global, paste0(path.diff, '/tree_ring_global_', topN, '.txt'), 
                    sep = '\t', quote = F, row.names = F, col.names = F, na = '')
        write.table(df.ring, paste0(path.diff, '/tree_ring_', topN, '.txt'), 
                    sep = '\t', quote = F, row.names = F, col.names = F, na = '')
        
    }
}


# shell 
# 0123
# female
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/0123/female_diff_1e-04
rm -rf track*
cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_female_0123.pdf --size 8

# male
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/0123/male_diff_1e-04
rm -rf track*
cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_male_0123.pdf --size 8

# 01
# female
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/01/female_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_female_01.pdf --size 8

# male
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/01/male_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_male_01.pdf --size 8

# 02
# female
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/02/female_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_female_02.pdf --size 8

# male
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/02/male_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_male_02.pdf --size 8


# 03
# female
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/03/female_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_female_03.pdf --size 8

# male
cd /home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan/03/male_diff_1e-04
rm -rf track*
    cat ../../cfg/global.cfg ../tree_label_color_50.txt > track0
cat tree_ring_global_50.txt tree_ring_50.txt > track1
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt ../tree_backbone_50.txt graphlan.xml
# 绘图，size决定图片大小，越大字越小
graphlan.py graphlan.xml graphlan50_tree_male_03.pdf --size 8


