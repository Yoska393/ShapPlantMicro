
##```{r}
pheno19 <- read.csv(here("2019", "shoot_blup_all_2019.csv"), row.names = 1)
plot19<-rownames(pheno19)
var19<-pheno19$var.id
##```

##```{r}
grm <- readRDS(here("data", "grm_2.rds"))
head(grm)[,1:3]
#```


#```{r}

bm<- readRDS(here("data/genome","genoMarker_LD0.9_SNP234781.RDS"))
bm2<- readRDS(here("data/genome","genoMarker_LD0.5_SNP64438.RDS"))
bm3<- readRDS(here("data/genome","genoMarker_LD0.1_SNP16419.RDS"))
bm4<- readRDS(here("data/genome","genoMarker_LD0.001_SNP3078.RDS"))

head(bm)[,1:3]

genome <- as.matrix(bm4)
dim(genome)
#```


#```{r}
met<- read.csv(here("2019","2019_Tottori_Jul_RootMetabolome_X.csv"))
rownames(met) <- paste0(substr(as.matrix(met[7]),1,2),".",sprintf("%03d",as.matrix(met[8])))
met19 <- met[,10:ncol(met)]
colnames(met19)<- sub("M","X",colnames(met19))

dim(met19)

#```


#```{r}
metind<- read.csv(here("data","Metabolome_index.csv"))
metindex<- matrix(metind[,2],583,1)
rownames(metindex)<- metind[1:583,1]

metind[metind == ""] <- NA
# Omit rows that contain only NA
metind <- metind[rowSums(is.na(metind)) < ncol(metind), ]

# Omit columns that contain only NA
metind <- metind[, colSums(is.na(metind)) < nrow(metind)]


Xname<-colnames(met19)
# Create a named vector from metind
annotation_lookup <- setNames(metind$Annotation, metind$Name)

# Replace column names in m2s based on annotation_lookup
metname <- annotation_lookup[colnames(met19)]
#colnames(met19.cd) <- metname
metname_vec <- unname(annotation_lookup[colnames(met19)])

metname<-colnames(met19)
new_colnames <- paste0("X", 1:ncol(met19))
colnames(met19) <- new_colnames


ID<-readRDS(here("2019", "ID_Taxonomy.RDS"))

com19<-readRDS(here("2019", "2019_Micro_Tung.RDS"))
colnames(com19)<- ID[,2]
rownames(com19)<- sub("-",".",rownames(com19))

com19<- com19[,!grepl("Mitochondria|Chloroplast", colnames(com19)) ]

comr <- read.csv(here("2019", "table-from-biom-7.csv"), row.names = 1)
comr  <- t(comr)

com19<- com19[rownames(com19) %in% rownames(comr),]
comr<- comr[rownames(comr) %in% rownames(com19),]


microname<-colnames(com19)
new_colnames <- paste0("B", 1:ncol(com19))
colnames(com19) <- new_colnames




### Drought Control

#```{r}
CD ="W4"
#```

#```{r}

if (CD == "W4"){
	print(CD)
	cd<- "Drought"
	
} else if (CD == "W1"){
	print(CD)
	cd<- "Control"
	
} else {
	print("error")
}

#```

#```{r}
selector.grm2 <- pheno19$var.id %in% rownames(grm)
selector.com2 <-rownames(pheno19) %in% rownames(com19)


pl<- list(pheno19)
pl.sel<- list()

for (i in 1){
	pheno<- pl[[i]]
	
	
	selector.grm <- selector.grm2
	selector.com <- selector.com2
	
	
	pheno.sel <- pheno[selector.grm & selector.com, ]
	#head(pheno.sel)
	print(dim(pheno.sel))
	
	selector.cd <- substr(rownames(pheno.sel), 1, 2) == CD 
	pheno.sel.cd <- pheno.sel[selector.cd, ]
	dim(pheno.sel.cd)
	
	plot.id <- rownames(pheno.sel.cd)
	rownames(pheno.sel.cd)<-pheno.sel.cd$var.id
	pheno.sel.cd$var.id <- plot.id
	colnames(pheno.sel.cd)[1]<- "plot.id"
	sg <- rownames(grm) %in% rownames(pheno.sel.cd)
	gv<-rownames(grm[sg,sg])
	pheno.sel.cd<- pheno.sel.cd[gv,]
	print(dim(pheno.sel.cd))
	pl.sel<- pheno.sel.cd
}

#psc18<-  pl.sel[[1]]
psc19<-  pl.sel

#s_p <- colnames(psc19) %in% colnames(psc20)
#psc19<- psc19[,s_p]
#s_p <- colnames(psc20) %in% colnames(psc19)
#psc20<- psc20[,s_p]

# table(rownames(psc19) == rownames(psc20))

#```


#```{r}
var.id <- rownames(psc19)
grm.cd <- grm[var.id, var.id]
plot.id  <- psc19$plot.id
com19.cd <- com19[plot.id,]
dim(com19.cd)
table(apply(com19.cd, 2, sum) > 0.5)
table(apply(com19.cd, 2, sd) > 0.5)


if (cd == "Drought"){
	CD<-"W4"
	com19.cd<-as.matrix(com19.cd[, apply(com19.cd, 2, sd) > 0]) 
} else if (cd=="Control"){
	CD<-"W1"
	com19.cd<-as.matrix(com19.cd[, apply(com19.cd, 2, sd) > 0.5]) 
}


# 9444 for W1 sd >0
#5424 for W1 (sd>0.5)
# 4946 for W4  (sd>0)

comr.cd <- comr[plot.id,]
dim(comr.cd)
table(apply(comr.cd, 2, sum) > 0)
table(apply(comr.cd, 2, sd) > 0)
comr.cd<-as.matrix(comr.cd[, apply(comr.cd, 2, sum) > 0]) 

#met19ns.cd <- met19[plot.id, ]
met19.cd <- met19[plot.id, ]


cm<-colSums(met19.cd)
cc<-cm>0.04
#table(cc)
met19.cd<-met19.cd[,cc ]


#met19.cd <- scale(met19[plot.id, ])
#predmet19ns<-predmet19
#predmet19<-scale(predmet19)
#predmet19<-predmet19
genome.cd<- genome[var.id,]
#predmet19RF<-scale(predmet19RF)
#predmet19f20RF<-scale(predmet19f20RF)

#predmet19RF<-predmet19RF






#```{r, message =F, include=TRUE,eval=TRUE}
m0  <- genome.cd #genome
m1  <- grm.cd #grm
m2  <- met19.cd #metabolome
m2s  <-scale(m2)
m3  <- com19.cd # microbiome 
m3s  <-scale(m3)
m4  <- comr.cd # microbiome 
m13 <- psc19[,-1] # phenotype


print(cd)
dim(m0)
dim(m1)
dim(m2)
dim(m3)
dim(m4)
dim(m13)

MMlist <- list(
	metabolome = m2s,
	microbiome = m3s,
	metname = metname_vec,  # This assumes 'metname' was derived from 'met19.cd'
	microname = microname # This is directly the column names from 'com19.cd'
)

MMlist2<-list(cd=cd,
							genome=m0,
							grm=m1,
							met=m2,
							mic=m3,
							mir=m4,
							pheno=m13,
							metname = metname_vec,  # This assumes 'metname' was derived from 'met19.cd'
							microname = microname) # This is directly the column names from 'com19.cd')

#

#saveRDS(MMlist2,here("data","SoyData_Control2.RDS"))
#saveRDS(MMlist2,here("data","SoyData_Drought2.RDS"))

