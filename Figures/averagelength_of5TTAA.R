


###########To extract chromosome info from gff file
gff.file <- "https://plasmodb.org/common/downloads/release-58/PknowlesiH/gff/data/PlasmoDB-58_PknowlesiH.gff"
header.lines <- readLines(gff.file, n = 100)
#The lines with the standard chromosomes start with "##sequence-region PKNH_" but is not followed by "archive" anywhere in the text.
#Select them. ^ represents the start of a line or string, When used as ^##, it matches the beginning of a line that starts with the characters "##"
ll <- header.lines[grepl(header.lines, pattern = "^##sequence-region PKNH_")]
#only include chromosomes rather than  API, MIT genome
#ll <- ll[57:70]
ll <- ll[1:70]
#split them by space, and create a data.frame
gg <- data.frame(do.call(rbind, strsplit(ll, split = " ")))
gg[,3] <- as.numeric(as.character(gg[,3]))
gg[,4] <- as.numeric(as.character(gg[,4]))

#######Total length of genome
sum(gg$X4)

cm <- read.xlsx("./Output/count_matrix/all/Pk_count_matrix_run1_2_3merged_essentialomeOnly.xlsx")
###To remove API and MIT
cm2 <- cm %>% dplyr::filter(!grepl("API", GeneID, fixed = TRUE) & !grepl("MIT", GeneID, fixed = TRUE))
dim(cm2)#319266

sum(gg$X4)/159633
sum(gg$X4)/159633*5
