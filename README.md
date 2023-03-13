## Course Assignments

A description of files within this folder:

-   `UNIX_Assignment.md` and `UNIX_Assignment.pdf`: Instructions for the assignment
-   `UNIX_Assignment_Template.md` and `UNIX_Assignment_Template.pdf`: An example of what your Markdown file should look like when you submit your assignment, including some Markdown syntax that should be helpful for you. The pdf shows how this file is rendered using a tool such as "MacDown"
-   The two files `fang_et_al_genotypes.txt` and `snp_positions.txt` are data files you will be reformatting for the assignment
-   The `transpose.awk` script will be needed to transpose the data (see instructions in `UNIX_Assignment.md`)

# BCB546-R-Assignememt

```{r}

library(tidyverse)
library(tidyr)
library(readr)

fang=read_tsv("fang_et_al_genotypes.txt", col_names=TRUE)
write_tsv(fang,"read_fang.tsv")

snp=read_tsv("snp_position.txt", col_names=TRUE)
write_tsv(snp,"read_snp_position.tsv")
```

# Data Inspection

```{r}

summary("read_snp_position.tsv")
summary("R_assignment/snp_position.txt")
summary("C:/Users/upatel/OneDrive - Iowa State University/BBMB Classes/BCB546/EEOB546_R_lesson/R_assignment/read_fang.tsv")
nrow("R_assignment/snp_position.txt")
nrow("R_assignment/read_snp_position.tsv")
ncol("R_assignment/read_snp_position.tsv")
dim("R_assignment/read_snp_position.tsv")
head("R_assignment/snp_position.txt")
ncol("R_assignment/read_fang.tsv")
view(read_tsv("read_snp_position.tsv"))
view(read_tsv("read_fang.tsv"))
```

# Data Processing

### Exrtraction of columns from SNP_position file

```{r}
threecolsnp=select(snp,1,3,4)
write_tsv(threecolsnp,"3columnssnp.tsv")
```

### Maize Processing

```{r}
data=read.delim("read_fang.tsv")
corn_fang= grepl("ZMMIL",data$Group)|grepl("ZMMLR",data$Group)|grepl("ZMMMR",data$Group)
result=data[corn_fang,]
write_tsv(result,"maize_fang.tsv")
```

### Transposing maize

```{r}
transp_maize= read.table("maize_fang.tsv", header = TRUE ,sep ="\t")
transposed_maize=t(transp_maize)
write.table(transposed_maize,"transposed_maize.tsv",sep= "\t",row.names=TRUE)
```

### Teosinte Processing

```{r}

data=read.delim("read_fang.tsv")
corn_fang= grepl("ZMPIL",data$Group)|grepl("ZMPBA",data$Group)|grepl("ZMPJA",data$Group)
result=data[corn_fang,]
write_tsv(result,"teosinte_fang.tsv")
```

### Transposing teosinte

```{r}
transp_teosinte= read.table("teosinte_fang.tsv", header = TRUE ,sep ="\t")
transposed_teosinte=t(transp_teosinte)
write.table(transposed_teosinte,"transposed_teosinte.tsv",sep= "\t",row.names=TRUE)
```

## Merging SNP with Maize and Teosinte

### Join SNP with Maize

```{r}
remove_header= read.table("transposed_maize.tsv", header= TRUE)
x=remove_header[-c(1:2),]
write.table(x,"transposed_maize1.tsv",sep= "\t")
maize= read.table("transposed_maize1.tsv",header=TRUE,sep="\t", row.names= 1)
snp3= read.table("3columnssnp.tsv", header=TRUE, sep="\t", row.names= 1)
snp_maize=merge(snp3,maize, by = 0)
write.table(snp_maize,"snp_maize.tsv", sep= "\t", row.names= TRUE)
```

### Join SNP with Teosinte

```{r}
remove_header= read.table("transposed_teosinte.tsv", header= TRUE)
x=remove_header[-c(1:2),]
write.table(x,"transposed_teosinte1.tsv",sep= "\t")
teosinte= read.table("transposed_teosinte1.tsv",header=TRUE,sep="\t", row.names= 1)
snp3= read.table("3columnssnp.tsv", header=TRUE, sep="\t", row.names= 1)
snp_teosinte=merge(snp3,teosinte, by = 0)
write.table(snp_teosinte,"snp_teosinte.tsv", sep= "\t", row.names= TRUE)
```

## Extracting Unknown and Multiple position in both Maize and Teosinte files

#### Maize Multiple and Unknown Extraction

```{r}
read= read.delim("snp_maize.tsv")
unkonwn_maize=grepl("unknown", read$Chromosome) | grepl("unknown", read$Position)
multiple_maize=grepl("multiple", read$Chromosome) | grepl("multiple",read$Position)
result_unknown=read[unkonwn_maize,]
result_multiple=read[multiple_maize,]
write_tsv(result_unknown,"unknown_maize.tsv")
write_tsv(result_multiple,"multiple_maize.tsv")
```

#### Teosinte Multiple and Unknown Extraction

```{r}
read= read.delim("snp_teosinte.tsv")
unkonwn_teosinte=grepl("unknown", read$Chromosome) | grepl("unknown", read$Position)
multiple_teosinte=grepl("multiple", read$Chromosome) | grepl("multiple",read$Position)
result_unknown=read[unkonwn_teosinte,]
result_multiple=read[multiple_teosinte,]
write_tsv(result_unknown,"unknown_teosinte.tsv")
write_tsv(result_multiple,"multiple_teosinte.tsv")
```

## Extraction of Chromosomes 1-10 for Maize and Teosinte

#### Maize Chromosomes 1-10 Extraction:

#### for writing a new file i have sub-directory for my reference for only Chromosome1, so it might through an error when you run the script.However, if you run from chromosome 2 it shoudl create file in you directory

```{r}
read_maize= read.delim("snp_maize.tsv")
Chr1_maize=grepl("^1$", read_maize$Chromosome)
reChr1=read[Chr1_maize,]
write_tsv(reChr1,"Chr1_maize.tsv")

Chr2_maize=grepl("^2$", read_maize$Chromosome)
reChr2=read[Chr2_maize,]
write_tsv(reChr2,"Chr2_maize.tsv")

Chr3_maize=grepl("^3$", read_maize$Chromosome)
reChr3=read[Chr3_maize,]
write_tsv(reChr3,"Chr3_maize.tsv")

Chr4_maize=grepl("^4$", read_maize$Chromosome)
reChr4=read[Chr4_maize,]
write_tsv(reChr4,"Chr4_maize.tsv")

Chr5_maize=grepl("^5$", read_maize$Chromosome)
reChr5=read[Chr5_maize,]
write_tsv(reChr5,"Chr5_maize.tsv")

Chr6_maize=grepl("^6$", read_maize$Chromosome)
reChr6=read[Chr6_maize,]
write_tsv(reChr6,"Chr6_maize.tsv")

Chr7_maize=grepl("^7$", read_maize$Chromosome)
reChr7=read[Chr7_maize,]
write_tsv(reChr7,"Chr7_maize.tsv")

Chr8_maize=grepl("^8$", read_maize$Chromosome)
reChr8=read[Chr8_maize,]
write_tsv(reChr8,"Chr8_maize.tsv")

Chr9_maize=grepl("^9$", read_maize$Chromosome)
reChr9=read[Chr9_maize,]
write_tsv(reChr9,"Chr9_maize.tsv")

Chr10_maize=grepl("^10$", read_maize$Chromosome)
reChr10=read[Chr10_maize,]
write_tsv(reChr10,"Chr10_maize.tsv")
```

#### Teosinte Chromosomes 1-10 Extraction

```{r}
read_teosinte= read.delim("snp_teosinte.tsv")
Chr1_teosinte=grepl("^1$", read_teosinte$Chromosome)
reChr1=read[Chr1_teosinte,]
write_tsv(reChr1,"Chr1_teosinte.tsv")

Chr2_teosinte=grepl("^2$", read_teosinte$Chromosome)
reChr2=read[Chr2_teosinte,]
write_tsv(reChr2,"Chr2_teosinte.tsv")

Chr3_teosinte=grepl("^3$", read_teosinte$Chromosome)
reChr3=read[Chr3_teosinte,]
write_tsv(reChr3,"Chr3_teosinte.tsv")

Chr4_teosinte=grepl("^4$", read_teosinte$Chromosome)
reChr4=read[Chr4_teosinte,]
write_tsv(reChr4,"Chr4_teosinte.tsv")

Chr5_teosinte=grepl("^5$", read_teosinte$Chromosome)
reChr5=read[Chr5_teosinte,]
write_tsv(reChr5,"Chr5_teosinte.tsv")

Chr6_teosinte=grepl("^6$", read_teosinte$Chromosome)
reChr6=read[Chr6_teosinte,]
write_tsv(reChr6,"Chr6_teosinte.tsv")

Chr7_teosinte=grepl("^7$", read_teosinte$Chromosome)
reChr7=read[Chr7_teosinte,]
write_tsv(reChr7,"Chr7_teosinte.tsv")

Chr8_teosinte=grepl("^8$", read_teosinte$Chromosome)
reChr8=read[Chr8_teosinte,]
write_tsv(reChr8,"Chr8_teosinte.tsv")

Chr9_teosinte=grepl("^9$", read_teosinte$Chromosome)
reChr9=read[Chr9_teosinte,]
write_tsv(reChr9,"Chr9_teosinte.tsv")

Chr10_teosinte=grepl("^10$", read_teosinte$Chromosome)
reChr10=read[Chr10_teosinte,]
write_tsv(reChr10,"Chr10_teosinte.tsv")
```

## Sorting based on Position and with "?" in increasing order for each chromosomes in Maize and Teosinte

#### for writing a new file i have sub-directory for my reference for only Chromosome1, so it might through an error when you run the script. However, if you run from chromosome 2 it shoudl create file in you directory

#### Maize

```{r}
read_maize= read.table("Chr1_maize.tsv", header = TRUE, sep = "\t")
maize_position1= arrange(read_maize,Position)
write.table(maize_position1, "incresing_maize_chr1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize2= read.table("Chr2_maize.tsv", header = TRUE, sep = "\t")
maize_position2= arrange(read_maize2,Position)
write.table(maize_position2, "incresing_maize_chr2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize3= read.table("Chr3_maize.tsv", header = TRUE, sep = "\t")
maize_position3= arrange(read_maize3,Position)
write.table(maize_position3, "incresing_maize_chr3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize4= read.table("Chr4_maize.tsv", header = TRUE, sep = "\t")
maize_position4= arrange(read_maize4,Position)
write.table(maize_position4, "incresing_maize_chr4.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize5= read.table("Chr5_maize.tsv", header = TRUE, sep = "\t")
maize_position5= arrange(read_maize5,Position)
write.table(maize_position5, "incresing_maize_chr5.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize6= read.table("Chr6_maize.tsv", header = TRUE, sep = "\t")
maize_position6= arrange(read_maize6,Position)
write.table(maize_position6, "incresing_maize_chr6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize7= read.table("Chr7_maize.tsv", header = TRUE, sep = "\t")
maize_position7= arrange(read_maize7,Position)
write.table(maize_position7, "incresing_maize_chr7.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize8= read.table("Chr8_maize.tsv", header = TRUE, sep = "\t")
maize_position8= arrange(read_maize8,Position)
write.table(maize_position8, "incresing_maize_chr8.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize9= read.table("Chr9_maize.tsv", header = TRUE, sep = "\t")
maize_position9= arrange(read_maize9,Position)
write.table(maize_position9, "incresing_maize_chr9.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_maize10= read.table("Chr10_maize.tsv", header = TRUE, sep = "\t")
maize_position10= arrange(read_maize10,Position)
write.table(maize_position10, "incresing_maize_chr10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

### Teosinte

```{r}
read_teosinte1= read.table("Chr1_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position1= arrange(read_teosinte1,Position)
write.table(teosinte_position1, "incresing_teosinte_chr1.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte2= read.table("Chr2_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position2= arrange(read_teosinte2,Position)
write.table(teosinte_position2, "incresing_teosinte_chr2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte3= read.table("Chr3_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position3= arrange(read_teosinte3,Position)
write.table(teosinte_position3, "incresing_teosinte_chr3.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte4= read.table("Chr4_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position4= arrange(read_teosinte4,Position)
write.table(teosinte_position4, "incresing_teosinte_chr4.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte5= read.table("Chr5_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position5= arrange(read_teosinte5,Position)
write.table(teosinte_position5, "incresing_teosinte_chr5.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte6= read.table("Chr6_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position6= arrange(read_teosinte6,Position)
write.table(teosinte_position6, "incresing_teosinte_chr6.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte7= read.table("Chr7_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position7= arrange(read_teosinte7,Position)
write.table(teosinte_position7, "incresing_teosinte_chr7.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte8= read.table("Chr8_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position8= arrange(read_teosinte8,Position)
write.table(teosinte_position8, "incresing_teosinte_chr8.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte9= read.table("Chr9_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position9= arrange(read_teosinte9,Position)
write.table(teosinte_position9, "incresing_teosinte_chr9.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

read_teosinte10= read.table("Chr10_teosinte.tsv", header = TRUE, sep = "\t")
teosinte_position10= arrange(read_teosinte10,Position)
write.table(teosinte_position10, "incresing_teosinte_chr10.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

## Sorting based on positon and "-" in decreading order for each chromosomes in Maize and Teosinte

#### for writing a new file i have sub-directory for my reference for only Chromosome1, so it might through an error when you run the script. However, if you run from chromosome 2 it shoudl create file in you directory.

#### Maize
```{r}
read_maize1= read.table("Chr1_maize.tsv", header = TRUE, sep = "\t")
maize1_position1= arrange(read_maize1,desc(Position))
maize1_position1[maize1_position1 == "?/?"] = "-/-"
write.table(maize1_position1, file = "decreasing_maize_chr1.tsv", sep = "\t",quote=FALSE,row.names = FALSE)

read_maize2= read.table("Chr2_maize.tsv", header = TRUE, sep = "\t")
maize2_position2= arrange(read_maize2,desc(Position))
maize2_position2[maize2_position2 == "?/?"] = "-/-"
write.table(maize2_position2, file = "decreasing_maize_chr2.tsv", sep = "\t", quote=FALSE, row.names = FALSE)

read_maize3= read.table("Chr3_maize.tsv", header = TRUE, sep = "\t")
maize3_position3= arrange(read_maize3,desc(Position))
maize3_position3[maize3_position3 == "?/?"] = "-/-"
write.table(maize3_position3, file = "decreasing_maize_chr3.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize4= read.table("Chr4_maize.tsv", header = TRUE, sep = "\t")
maize4_position4= arrange(read_maize4,desc(Position))
maize4_position4[maize4_position4 == "?/?"] = "-/-"
write.table(maize4_position4, file = "decreasing_maize_chr4.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize5= read.table("Chr5_maize.tsv", header = TRUE, sep = "\t")
maize5_position5= arrange(read_maize5,desc(Position))
maize5_position5[maize5_position5 == "?/?"] = "-/-"
write.table(maize5_position5, file = "decreasing_maize_chr5.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize6= read.table("Chr6_maize.tsv", header = TRUE, sep = "\t")
maize6_position6= arrange(read_maize6,desc(Position))
maize6_position6[maize6_position6 == "?/?"] = "-/-"
write.table(maize6_position6, file = "decreasing_maize_chr6.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize7= read.table("Chr7_maize.tsv", header = TRUE, sep = "\t")
maize7_position7= arrange(read_maize7,desc(Position))
maize7_position7[maize7_position7 == "?/?"] = "-/-"
write.table(maize7_position7, file = "decreasing_maize_chr7.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize8= read.table("Chr8_maize.tsv", header = TRUE, sep = "\t")
maize8_position8= arrange(read_maize8,desc(Position))
maize8_position8[maize8_position8 == "?/?"] = "-/-"
write.table(maize8_position8, file = "decreasing_maize_chr8.tsv", sep = "\t",quote=FALSE, row.names = FALSE)


read_maize9= read.table("Chr9_maize.tsv", header = TRUE, sep = "\t")
maize9_position9= arrange(read_maize9,desc(Position))
maize9_position9[maize9_position9 == "?/?"] = "-/-"
write.table(maize9_position9, file = "decreasing_maize_chr9.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_maize10= read.table("Chr10_maize.tsv", header = TRUE, sep = "\t")
maize10_position10= arrange(read_maize10,desc(Position))
maize10_position10[maize10_position10 == "?/?"] = "-/-"
write.table(maize10_position10, file = "decreasing_maize_chr10.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

```

### Teosinte
```{r}
read_teosinte1= read.table("teosintechr/Chr1_teosinte.tsv", header = TRUE, sep = "\t")
teosinte1_position1= arrange(read_teosinte1,desc(Position))
teosinte1_position1[teosinte1_position1 == "?/?"] = "-/-"
write.table(teosinte1_position1, file = "teosintechr/decreasing_teosinte_chr1.tsv", sep = "\t",quote=FALSE,row.names = FALSE)

read_teosinte2= read.table("Chr2_teosinte.tsv", header = TRUE, sep = "\t")
teosinte2_position2= arrange(read_teosinte2,desc(Position))
teosinte2_position2[teosinte2_position2 == "?/?"] = "-/-"
write.table(teosinte2_position2, file = "decreasing_teosinte_chr2.tsv", sep = "\t", quote=FALSE, row.names = FALSE)

read_teosinte3= read.table("Chr3_teosinte.tsv", header = TRUE, sep = "\t")
teosinte3_position3= arrange(read_teosinte3,desc(Position))
teosinte3_position3[teosinte3_position3 == "?/?"] = "-/-"
write.table(teosinte3_position3, file = "decreasing_teosinte_chr3.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte4= read.table("Chr4_teosinte.tsv", header = TRUE, sep = "\t")
teosinte4_position4= arrange(read_teosinte4,desc(Position))
teosinte4_position4[teosinte4_position4 == "?/?"] = "-/-"
write.table(teosinte4_position4, file = "decreasing_teosinte_chr4.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte5= read.table("Chr5_teosinte.tsv", header = TRUE, sep = "\t")
teosinte5_position5= arrange(read_teosinte5,desc(Position))
teosinte5_position5[teosinte5_position5 == "?/?"] = "-/-"
write.table(teosinte5_position5, file = "decreasing_teosinte_chr5.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte6= read.table("Chr6_teosinte.tsv", header = TRUE, sep = "\t")
teosinte6_position6= arrange(read_teosinte6,desc(Position))
teosinte6_position6[teosinte6_position6 == "?/?"] = "-/-"
write.table(teosinte6_position6, file = "decreasing_teosinte_chr6.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte7= read.table("Chr7_teosinte.tsv", header = TRUE, sep = "\t")
teosinte7_position7= arrange(read_teosinte7,desc(Position))
teosinte7_position7[teosinte7_position7 == "?/?"] = "-/-"
write.table(teosinte7_position7, file = "decreasing_teosinte_chr7.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte8= read.table("Chr8_teosinte.tsv", header = TRUE, sep = "\t")
teosinte8_position8= arrange(read_teosinte8,desc(Position))
teosinte8_position8[teosinte8_position8 == "?/?"] = "-/-"
write.table(teosinte8_position8, file = "decreasing_teosinte_chr8.tsv", sep = "\t",quote=FALSE, row.names = FALSE)


read_teosinte9= read.table("Chr9_teosinte.tsv", header = TRUE, sep = "\t")
teosinte9_position9= arrange(read_teosinte9,desc(Position))
teosinte9_position9[teosinte9_position9 == "?/?"] = "-/-"
write.table(teosinte9_position9, file = "decreasing_teosinte_chr9.tsv", sep = "\t",quote=FALSE, row.names = FALSE)

read_teosinte10= read.table("Chr10_teosinte.tsv", header = TRUE, sep = "\t")
teosinte10_position10= arrange(read_teosinte10,desc(Position))
teosinte10_position10[teosinte10_position10 == "?/?"] = "-/-"
write.table(teosinte10_position10, file = "decreasing_teosinte_chr10.tsv", sep = "\t",quote=FALSE, row.names = FALSE)
```

# Visulaization

```{r}
install.packages("ggplot2")

```

```{r}
```

```{r}
```
