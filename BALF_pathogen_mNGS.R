
###ABFV loading and pre-analysis
ABFV <- read_tsv("Downloads/taxonomy.bracken.ABFV.txt",num_threads = 20) %>%
        filter(rank=="S") %>% filter(!(taxid %in% c(77133,2832643))) 
colnames(ABFV)=gsub(".bracken.new","",colnames(ABFV))
ABFV_s <- ABFV  %>%filter(rank=="S")   %>%filter(!(taxid %in% c(77133,2832643))) %>% 
        select(-rank,-taxid) %>% column_to_rownames("tax_name")
#ABFV_s <- ABFV %>%filter(rank=="S")%>% filter(taxid %in% (filter(abfv_list,X1=="Fungi"|X1=="Bacterial") %>% pull(X2))) %>%select(-rank,-taxid) %>% column_to_rownames("tax_name")

nc_s <- read_tsv("Downloads/Blood.bracken.ABFV.txt",num_threads = 5) %>% filter(rank=="S") %>% filter(!(taxid %in% c(77133,2832643))) %>% select(-rank,-taxid) %>%
        column_to_rownames("tax_name") %>% as.matrix() %>% prop.table(.,2) %>%.[rowSums(.>0.005)>0,] %>% prop.table(.,2)
nc_meta <- read_tsv("Downloads/nc_metadata.tsv",)
nc_meta$date=as.Date(nc_meta$date, format = "%m/%d/%Y")
nc_meta$group <- str_sub(nc_meta$id, start = -1,end = -1)
find_closest <- function(row, matrix2) {
  group_match <- matrix2 %>% filter(Strategy == row$Strategy)
  if (nrow(group_match) == 0) return(data.frame(ID = NA))
  closest_row <- group_match %>%
    mutate(date_diff = abs(difftime(date, row$Date, units = "days"))) %>%
    arrange(date_diff) %>%
    slice(1) %>% pull(id)
  return(closest_row)
}
cc=list()
tt=data.frame(ID=colnames(ABFV_s),counts=colSums(ABFV_s)) %>% mutate(group=str_sub(ID,-1),UID=gsub("-.*","",ID))
for (i in tt$ID) {
  method=str_sub(i,-1)
  df=ABFV_s_p[,i,drop=F] %>% as.data.frame() %>% rownames_to_column("tax") %>%
    full_join(.,nc_s %>% as.data.frame()%>% dplyr::select(ends_with(method)) %>%
                rownames_to_column("tax"),)%>% as.data.frame %>%  column_to_rownames("tax")
  df[is.na(df)]=0
  df1=JSD(t(df))
  colnames(df1) <- colnames(df)
  rownames(df1) <- colnames(df)
  tt[i,"nc"]=df1[,i] %>% as.data.frame() %>% filter(`.`!=0) %>% arrange(desc(`.`)) %>% tail(1) %>% rownames()
}
metadata$nc= apply(metadata,1,function(x){find_closest(as.data.frame(t(x)),nc_meta)})
metadata$nc <-  as.character(metadata$nc)

###Decontam based on 5-fold nc
tt1=filter(tt,ID %in% ll_pair)
tt1$pair=NA
for (id in unique(tt1$UID)) {
  ll1 <- str_subset(ll_pair,paste0(id,"-"))
  tt1[ll1[1],"pair"]=ll1[2]
  tt1[ll1[2],"pair"]=ll1[1]
}


ll1 <- str_subset(ll_pair,paste0(id,"-"))
lr <- ll1[endsWith(ll1,"R")]
ld <- ll1[endsWith(ll1,"D")]


ll_spe=rownames(ABFV_s_p)
ll=filter(tt1,nc==i) %>% pull(ID)
df1=ABFV_s_p[,ll]
ll_except=nc_s[,i ,drop=F] %>% filter(get(i)>0) %>% rownames(.) %>% intersect(ll_spe,.)
for (k in ll_except) {
  j=5*nc_s[k,i]
  df1[k,][df1[k,]<j]=0
}
dd=df1[ll_except,] %>% rownames_to_column("tax")%>% melt() %>% filter(value!=0)
dd$variable=as.character(dd$variable)
for (i in 1:nrow(dd) ) {
  if(ABFV_s_p[dd[i,"tax"],tt1[dd[i,"variable"],"pair"]]==0) {
    df1[dd[i,"tax"],dd[i,"variable"]]=0
  }
}
cc=list()
for (jj in unique(tt1$nc)) {
  ll=filter(tt1,nc==jj) %>% pull(ID)
  df1=ABFV_s_p[,ll]
  ll_except=nc_s[,jj ,drop=F] %>% filter(get(jj)>0) %>% rownames(.) %>% intersect(ll_spe,.)
  for (k in ll_except) {
    j=5*nc_s[k,jj]
    df1[k,][df1[k,]<j]=0
  }
  dd=df1[ll_except,] %>% rownames_to_column("tax")%>% melt() %>% filter(value!=0)
  dd$variable=as.character(dd$variable)
  for (i in 1:nrow(dd) ) {
    if(ABFV_s_p[dd[i,"tax"],tt1[dd[i,"variable"],"pair"]]==0) {
      df1[dd[i,"tax"],dd[i,"variable"]]=0
    }
  }
cc[[jj]]=df1
  }
df=do.call(cbind,cc)
colnames(df)=gsub(".*\\.","",colnames(df))
df1=ABFV_s[,colnames(df)]
df1=df1[rownames(df),]
df1[df==0]=0
tt1$qc <-  df1[,tt1$ID] %>% colSums()
ABFV_s_de=df1



add_meta <- add_meta %>% filter(!is.na(ID))
rownames(add_meta)=add_meta$ID
metadata <- rbind(metadata,add_meta)
metadata <- metadata %>%
        mutate(
                Age_group = case_when(
                        Age <= 0.5 ~ "Infant group",
                        Age > 0.5 & Age <= 2 ~ "Toddler group",
                        #Age <= 2 ~ "Infant group",
                        Age > 2 & Age <= 6 ~ "Children group",
                        Age > 6 & Age <= 18 ~ "Teenager group",
                        Age > 18 & Age <= 35 ~ "Youth group",
                        Age > 35 & Age <= 60 ~ "Middle-aged group",
                        Age > 60 ~ "Senior group",
                        TRUE ~ NA_character_  # 
                        )) %>% mutate(Age_group=factor(Age_group,levels=c("Infant group","Toddler group","Children group","Teenager group","Youth group","Middle-aged group","Senior group")))
metadata$month <-  gsub("/.*","",metadata$Date) %>% as.numeric()
metadata <- metadata %>%
        mutate(
                season = case_when(
                        month %in% c(3, 4, 5) ~ "Spring",
                        month %in% c(6, 7, 8) ~ "Summer",
                        month %in% c(9, 10, 11) ~ "Autumn",
                        month %in% c(12, 1, 2) ~ "Winter",
                        TRUE ~ NA_character_  # 
                )) %>% mutate(season=factor(season,levels=c("Spring","Summer","Autumn","Winter")))
metadata$Sexual <- gsub("女","female",metadata$Sexual)
metadata$Sexual <- gsub("男","male",metadata$Sexual)
metadata$year <- gsub(".*/","",metadata$Date)
north_provinces <- c("北京", "天津", "河北", "山西", "内蒙古", "辽宁", "吉林", "黑龙江", "山东", "陕西", "宁夏", "甘肃", "新疆","河南","青海")
south_provinces <- c("上海", "江苏", "浙江", "安徽", "福建", "江西", "湖北", "湖南", "广东", "广西", "海南", "重庆", "四川", "贵州", "云南", "西藏")
metadata$south_north <- ifelse(metadata$Province %in% north_provinces, "Northern", 
                                ifelse(metadata$Province %in% south_provinces, "Southern", "未知"))

library(gtsummary)
filter(metadata,!duplicated(UID), UID %in% (colnames(ABFV) %>% gsub("-.*","",.) %>% unique() )) %>% 
        select(Sexual,Age,CRP  ,PCT,   WBC)%>%
        tbl_summary(digits = all_continuous() ~ 3)%>% as_flex_table() 





### RBM identify pathogen
find_species_to_split <- function(sample) {
        
        # species rank based on abundance
        sorted_sample <- sample[order(sample, decreasing = TRUE)]
        
        # find the biggest gap
        diff_max <- which.max(-diff(sorted_sample))
        
        names(sorted_sample)[1:diff_max]
}

find_rbm_spe <- function(sample) {
        ll1 <- str_subset(ll_pair,paste0(sample,"-"))
        if (all(ll1 %in% colnames(abfv_ss))) {
                dd=abfv_ss[order(rowSums(abfv_ss[,ll1]),decreasing = T),ll1] %>% head(100) %>% dist() %>% as.matrix()
                n <- nrow(dd)
                sorted_sample <- sapply(1:(n-1), function(i) dd[i+1, i])
                # 找到丰度差值最大的地方
                diff_max <- which.max(sorted_sample)
                print(sample)
                # 返回分隔点之前的物种
                rownames(dd)[1:diff_max]
        }

}

# apply RBM to every paired samples
ll=filter(abfv_list1,ABFV=="Virus") %>% pull(spe)
abfv_ss <- ABFV_s_p_de[setdiff(rownames(ABFV_s_p_de),ll),] %>% as.matrix() %>% prop.table(.,2) %>% as.data.frame()
ll=gsub("-.*","",ll_pair) %>% unique()
species_b <- sapply(ll,find_rbm_spe)
ll=filter(abfv_list1,ABFV=="Virus") %>% pull(spe)
abfv_ss <- ABFV_s_p_de[intersect(rownames(ABFV_s_p_de),ll),] %>% as.matrix() %>% prop.table(.,2)  %>% as.data.frame()
species_v_d <- apply(abfv_ss %>% select(ends_with("D")), 2, find_species_to_split)
species_v_r <- apply(abfv_ss %>% select(ends_with("R")), 2, find_species_to_split)

ll1 <- intersect(colnames(ABFV_s_p),ll_pair) %>% gsub("-.*","",.) %>% 
        table %>% as.data.frame() %>% filter(Freq==2) %>% pull(1) %>% as.character()
ll <- intersect(colnames(ABFV_s_p),ll_pair)%>% as.data.frame() %>% mutate(id=gsub("-.*","",.)) %>% filter(id %in% ll1) %>% pull(1)
df=select(ABFV_s_p[,ll],ends_with("D"))+select(ABFV_s_p[,ll],ends_with("R"))
species_f <- apply(df, 2, find_species_to_split)

pathogen <- read_tsv("Downloads/all_pathogen1.stats.tsv",col_names = F)
add_p <- read_tsv("Downloads/add.pathogen.tsv",col_names = F)
colnames(pathogen) <- c("ID","spe","kraken","genome","reads","YN","gap")
pathogen$ID=gsub("/pathogene.stat.tsv","",pathogen$ID)
pathogen$spe=gsub("_"," ",pathogen$spe)





ll=c("Acinetobacter baumannii","Pseudomonas aeruginosa","Orthopneumovirus hominis","Klebsiella pneumoniae","Corynebacterium striatum","Pneumocystis jirovecii","Candida albicans","Human alphaherpesvirus 1","Rhinovirus A","Haemophilus influenzae","Mycoplasmoides pneumoniae","Respirovirus pneumoniae","Staphylococcus aureus","Streptococcus pseudopneumoniae","Metapneumovirus hominis","Rhinovirus C","Stenotrophomonas maltophilia","Ralstonia pickettii","Mycobacterium tuberculosis","Streptococcus pneumoniae","Candida tropicalis","Chlamydia psittaci","Prevotella melaninogenica","Rothia mucilaginosa","Human betaherpesvirus 5","Enterococcus faecium","Moraxella catarrhalis","Escherichia coli","Tropheryma whipplei","Rhinovirus B","Influenza B virus","Staphylococcus epidermidis","Aspergillus fumigatus","Human coronavirus 229E","Respirovirus laryngotracheitidis","Pseudomonas fluorescens","Legionella pneumophila","Staphylococcus haemolyticus","Weissella viridescens","Chlamydia trachomatis","Betacoronavirus 1","Veillonella parvula","Bradyrhizobium sp. PSBB068")
dd[ll,] %>%as.data.frame() %>%  rownames_to_column("tax") %>%
        melt %>% na.omit() %>% mutate(value=log10(value+1e-5),group=str_sub(variable,-1)) %>% 
        gghistogram(.,"value",fill = "group",
                    palette = c("#00AFBB", "#E7B800"), rug = TRUE,color = "group")+
        facet_wrap(tax~.)+ 
        scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                      labels = trans_format("log10", math_format(10^.x)))


ll_pathogen=c("Human coronavirus NL63","Human coronavirus HKU1","Acinetobacter baumannii","Aspergillus flavus","Aspergillus fumigatus","Aspergillus niger","Aspergillus terreus","Bacteroides fragilis","Blastomyces dermatitidis","Bordetella pertussis","Burkholderia cepacia","Burkholderia pseudomallei","Chlamydia pneumoniae","Chlamydia psittaci","Citrobacter freundii","Citrobacter koseri","Coccidioides immitis","Coccidioides posadasii","Coxiella burnetii","Cryptococcus gattii VGI","Cryptococcus neoformans","Enterobacter cloacae","Escherichia coli","Francisella tularensis","Fusobacterium necrophorum","Fusobacterium nucleatum","Haemophilus influenzae","Histoplasma capsulatum","Human mastadenovirus B","Severe acute respiratory syndrome-related coronavirus","Metapneumovirus hominis","Human orthorubulavirus 4","Rhinovirus A","Influenza A virus","Klebsiella aerogenes","Klebsiella oxytoca","Klebsiella pneumoniae","Legionella pneumophila","Moraxella catarrhalis","Morganella morganii","Mycobacterium tuberculosis","Mycoplasmoides pneumoniae","Nocardia asteroides","Pasteurella multocida","Pneumocystis jirovecii","Proteus mirabilis","Pseudomonas aeruginosa","Orthopneumovirus hominis","Serratia marcescens","Staphylococcus aureus","Stenotrophomonas maltophilia","Streptococcus pneumoniae","Streptococcus pyogenes","Candida albicans","Respirovirus pneumoniae","Rhinovirus C","Respirovirus laryngotracheitidis","Rhinovirus B","Human coronavirus 229E","Human alphaherpesvirus 1","Human betaherpesvirus 5","Human gammaherpesvirus 4","Human mastadenovirus C","Influenza B virus","Influenza C virus","Enterococcus faecium","Betacoronavirus 1","Mucor piriformis","Corynebacterium striatum","Candida tropicalis","Tropheryma whipplei","Chlamydia trachomatis","Staphylococcus haemolyticus","Weissella viridescens")

ll1 <- dd2 %>% as.data.frame()%>%
        rownames_to_column("taxname") %>% melt %>% mutate(group=str_sub(variable,-1)) %>% group_by(taxname,group) %>% summarise(value=max(value)) %>% 
        arrange(desc(value)) %>% as.data.frame() %>% 
        filter(value>0.01) %>% pull(taxname) %>%
        unique() %>% intersect(.,ll_pathogen)
ll_nc_spe <- dd2 %>% as.data.frame()%>%
        rownames_to_column("taxname") %>% melt %>% mutate(group=str_sub(variable,-1)) %>% group_by(taxname,group) %>% summarise(value=max(value)) %>% 
        arrange(desc(value)) %>% as.data.frame() %>% 
        filter(value>0.01) %>% pull(taxname) %>%
        unique() %>% setdiff(.,ll_pathogen) %>% intersect(.,rownames(ABFV_s))
ll_nc_path <- dd2 %>% as.data.frame()%>%
        rownames_to_column("taxname") %>% melt %>% mutate(group=str_sub(variable,-1)) %>% group_by(taxname,group) %>% summarise(value=max(value)) %>% 
        arrange(desc(value)) %>% as.data.frame() %>% 
        filter(value>0.01) %>% pull(taxname) %>%
        unique() %>% intersect(.,ll_pathogen)


ll=colSums(ABFV_s_de) %>% as.data.frame() %>% .[colnames(ABFV_s_p),,drop=F] %>% rename("reads"=".") %>%
        rownames_to_column("ID")%>% mutate(group=str_sub(ID,-1),
                                           UID=gsub("-.*","",ID)) %>%
        group_by(UID) %>% mutate(un=length(unique(group)),nn=n()) %>%
        filter(un!=1) %>%filter(nn>2) %>%
        group_by(UID,group) %>% mutate(ii=n(),ui=paste(UID,group,sep="_")) %>%
  filter(ii>1) %>% arrange(ui,desc(reads)) %>% 
        filter(duplicated(ui)) %>% pull(ID)

ll1 <- colSums(ABFV_s_de) %>% as.data.frame() %>% .[colnames(ABFV_s_p),,drop=F] %>% rename("reads"=".") %>%
        rownames_to_column("ID")%>% mutate(group=str_sub(ID,-1),
                                           UID=gsub("-.*","",ID)) %>%
        group_by(UID) %>% mutate(un=length(unique(group)),nn=n()) %>%
        filter(un!=1) %>%filter(nn>2) %>% filter(!(ID %in% ll)) %>% pull(ID)



ll_pair=colSums(ABFV_s_de) %>% as.data.frame() %>% .[colnames(ABFV_s_p),,drop=F] %>% rename("reads"=".") %>%
        rownames_to_column("ID")%>% mutate(group=str_sub(ID,-1),
                                           UID=gsub("-.*","",ID)) %>%
        group_by(UID) %>% mutate(un=length(unique(group)),nn=n()) %>%
        filter(un!=1) %>%filter(nn==2) %>% pull(ID) %>% c(.,ll1)





#### pathogen identify function
df1 <- species_b %>% do.call(c,.) %>% as.data.frame() %>%
        rownames_to_column("ID") %>% mutate(ID=gsub("[0-9]$","",ID)) %>%
        rename("spe"=".") %>% right_join(df,.) 

ll=gsub("-.*","",ll_pair) %>% unique()
column_names <- c("UID", "spe", "krakenD", "krakenR", "gapD", "gapR", "bowtieD", "bowtieR", "RBM","Virus")

# 创建一个空的DataFrame
df <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(df) <- column_names
i=1
extract_info <- function(id, subset) {
        if (id %in% names(subset)) {
                return(subset[[id]])
        } else {
                return(NA)
        }
}
extract_gap <- function(id, subset) {
        if (id %in% names(subset)) {
                return(subset[[id]])
        } else {
                return(NA)
        }
}
i=1
for (id in ll) {
        ll1 <- str_subset(ll_pair,paste0(id,"-"))
        lr <- ll1[endsWith(ll1,"R")]
        ld <- ll1[endsWith(ll1,"D")]
        l1 <- unique(c(extract_info(ld,species_v_d),extract_info(lr,species_v_r),extract_info(id,species_b)))
        ll2 <- unique(c(filter(pathogen,ID %in% ll1) %>% pull(spe),l1)) %>% na.omit() %>% as.character()
        for (ii in ll2) {
          df[i,"UID"] = id
          df[i,"spe"] = ii
          df[i,"krakenD"] = ABFV_s_de[ii,ld]
          df[i,"krakenR"] = ABFV_s_de[ii,lr]
          if ((filter(pathogen,ID==ld,spe==ii) %>% nrow())>0) {
                  df[i,"krakenD"] = filter(pathogen,ID==ld,spe==ii) %>% pull(kraken)
          }
          if ((filter(pathogen,ID==lr,spe==ii) %>% nrow())>0) {
                  df[i,"krakenR"] = filter(pathogen,ID==lr,spe==ii) %>% pull(kraken)
          }
          df[i,"gapD"] = ifelse((filter(pathogen,ID==ld,spe==ii) %>% nrow())==0,NA,filter(pathogen,ID==ld,spe==ii) %>%pull(gap))
          df[i,"gapR"] = ifelse((filter(pathogen,ID==lr,spe==ii) %>% nrow())==0,NA,filter(pathogen,ID==lr,spe==ii) %>%pull(gap))
          df[i,"bowtieD"] = ifelse((filter(pathogen,ID==ld,spe==ii) %>% nrow())==0,NA,filter(pathogen,ID==ld,spe==ii) %>%pull(reads))
          df[i,"bowtieR"] = ifelse((filter(pathogen,ID==lr,spe==ii) %>% nrow())==0,NA,filter(pathogen,ID==lr,spe==ii) %>%pull(reads))
          df[i,"RBM"] =   ifelse((ii %in% l1)|grepl("virus",ii),"Y","N" )
          df[i,"Virus"] = ifelse(grepl("virus",ii),"Y","N")
          i=i+1
        }
        print(id)
}
#RBM <- df

df=filter(RBM,krakenD!=0|krakenR!=0) 
ll=c(species_v_d %>% do.call(c,.),species_v_r %>% do.call(c,.)) %>% table() %>% as.data.frame() %>% arrange(Freq) %>% filter(!grepl("virus",`.`)) 
df <- filter(df, !(spe %in% ll$.))
df <- filter(df, !(spe %in% ll))
df1 <- filter(df,Virus=="Y"|(krakenD!=0&krakenR!=0)) 
df1 <- filter(df1,Virus=="Y"|(krakenD!=0&krakenR!=0)) %>% filter((is.na(bowtieD)&is.na(bowtieR))|(Virus=="Y"&(gapD>0|gapR>0))|(gapD>1|gapR>1))
ll_pathogen2 <- c("Acinetobacter baumannii","Aspergillus flavus","Aspergillus fumigatus","Aspergillus niger","Aspergillus terreus","Bacteroides fragilis","Blastomyces dermatitidis","Bordetella pertussis","Burkholderia cepacia","Burkholderia pseudomallei","Chlamydia pneumoniae","Chlamydia psittaci","Citrobacter freundii","Citrobacter koseri","Coccidioides immitis","Coccidioides posadasii","Coxiella burnetii","Cryptococcus gattii VGI","Cryptococcus neoformans","Enterobacter cloacae","Escherichia coli","Francisella tularensis","Fusobacterium necrophorum","Fusobacterium nucleatum","Haemophilus influenzae","Histoplasma capsulatum","Human mastadenovirus B","Severe acute respiratory syndrome-related coronavirus","Metapneumovirus hominis","Human orthorubulavirus 4","Rhinovirus A","Influenza A virus","Klebsiella aerogenes","Klebsiella oxytoca","Klebsiella pneumoniae","Legionella pneumophila","Moraxella catarrhalis","Morganella morganii","Mycobacterium tuberculosis","Mycoplasmoides pneumoniae","Nocardia asteroides","Pasteurella multocida","Pneumocystis jirovecii","Proteus mirabilis","Pseudomonas aeruginosa","Orthopneumovirus hominis","Serratia marcescens","Staphylococcus aureus","Stenotrophomonas maltophilia","Streptococcus pneumoniae","Streptococcus pyogenes","Candida albicans","Respirovirus pneumoniae","Rhinovirus C","Respirovirus laryngotracheitidis","Rhinovirus B","Human coronavirus 229E","Human alphaherpesvirus 1","Human betaherpesvirus 5","Human gammaherpesvirus 4","Metapneumovirus hominis","Human mastadenovirus C","Influenza B virus","Influenza C virus","Enterococcus faecium","Betacoronavirus 1","Mucor piriformis","Candida tropicalis","Chlamydia trachomatis","Human coronavirus HKU1","Human betaherpesvirus 7","Human betaherpesvirus 6A","Human alphaherpesvirus 3","Enterovirus D","Human coronavirus NL63")
df2 <- filter(df1,RBM=="Y"|Virus=="Y") %>%
        mutate(pathogen= ifelse(spe %in% ll_pathogen2,"Y","N"))
df3=df2%>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD!=0&krakenR!=0))|
                          ((is.na(bowtieD)&is.na(bowtieR))|(gapD>1|gapR>1))|
                          grepl("Rhinovirus",spe)) 
df3u=anti_join(df1,df3)
df4 <- rbind(df3,df3u %>% filter(grepl("Mycobacterium tuberculosis",spe)) %>% filter(krakenD>5&krakenR>5) %>% mutate(pathogen= ifelse(spe %in% ll_pathogen2,"Y","N")))
df4    <-  right_join(metadata %>% filter(!duplicated(UID)),df4) %>% mutate(type=NA)
df4 <- right_join(abfv_list1,df4) %>% as.data.frame()
dd=filter(df4,is.na(gapD),is.na(gapR))
ll=table(dd$spe) %>% sort %>% as.data.frame %>% filter(grepl("virus",Var1))  %>% pull(Var1) %>% as.character() %>% setdiff(.,ll_pathogen2)
ll=ll[-which(ll=="Human betaherpesvirus 6B")]
df4 <- filter(df4,!(spe %in% ll))
ll=unique(df4$UID)

df1u <- filter(df1,!(UID %in% ll)) %>% filter(spe %in% ll_pathogen2) %>%
  filter(((is.na(bowtieD)&is.na(bowtieR))|(gapD>1|gapR>1))) %>% filter(krakenD>100,krakenR>100) %>% 
        filter(spe %in% c("Legionella pneumophila")) %>% mutate(pathogen= ifelse(spe %in% ll_pathogen2,"Y","N"))
df1u$RBM="Y"
df3=rbind(df1u,df2)%>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>5&krakenR>5))|
                          ((is.na(bowtieD)&is.na(bowtieR))|(gapD>1|gapR>1))|
                          grepl("Rhinovirus",spe)) 
df3u=anti_join(df1,df3)
df4 <- rbind(df3,df3u %>% filter(grepl("Mycobacterium tuberculosis",spe)) %>% filter(krakenD>3&krakenR>3) %>%
               mutate(pathogen= ifelse(spe %in% ll_pathogen2,"Y","N"))) 
ll1=filter(RBM,!(UID %in% unique(df4$UID))) %>% pull(UID) %>% unique()
ll2=setdiff(gsub("-.*","",ll_pair) %>% unique(),unique(df4$UID))
df4 <-        right_join(metadata %>% filter(!duplicated(UID)),df4) %>% mutate(type=NA)
df4 <- right_join(abfv_list1,df4) %>% as.data.frame()
dd=filter(df4,is.na(gapD),is.na(gapR))
ll=table(dd$spe) %>% sort %>% as.data.frame %>% filter(grepl("virus",Var1))  %>% pull(Var1) %>% as.character() %>% setdiff(.,ll_pathogen2)
df4 <- filter(df4,!(spe %in% ll))
df4$RBM="Y"
ll1=filter(RBM,!(UID %in% unique(df4$UID))) %>% pull(UID) %>% unique()


df4u <- right_join(metadata %>% filter(!duplicated(UID)),filter(RBM,UID %in% ll1) %>%
                            filter(RBM=="Y",Virus=="N") %>% filter(!duplicated(UID)) %>% 
                            mutate(pathogen= ifelse(spe %in% ll_pathogen2,"Y","N")) ) %>% mutate(type=NA)
df4u <- right_join(abfv_list1,df4u) %>% as.data.frame()
df4=rbind(df4,df4u)

### dedup Corona and rh  flu
df=filter(df4,grepl("corona",spe)) %>% filter(pathogen=="Y") %>% select(spe,UID,krakenD,krakenR,gapR,bowtieR,type)
ll=df$UID[duplicated(df$UID)]
for (i in ll) {
        l1 <- filter(df,UID==i) %>% arrange(krakenR) %>% pull(spe) %>% .[1]
        df4[which(df4$UID==i & df4$spe==l1),"pathogen"]="N"
}
df=filter(df4,grepl("Rhinovirus",spe)) %>% filter(pathogen=="Y") %>% select(spe,UID,krakenD,krakenR,gapR,bowtieR,type)
ll=df$UID[duplicated(df$UID)]
for (i in ll) {
        l1 <- filter(df,UID==i) %>% arrange(krakenR) %>% pull(spe) %>% .[1]
        df4[which(df4$UID==i & df4$spe==l1),"pathogen"]="N"
}

df4[grepl("herpesvirus",df4$spe),"pathogen"]="N"
df=filter(df4,grepl("Candida",spe)) %>% filter(pathogen=="Y") %>% select(spe,UID,krakenD,krakenR,gapR,bowtieR,type)

## exclude HHV
df4$spe <- gsub("Human gammaherpesvirus 4","EBV",df4$spe)
df4$spe <- gsub("Human betaherpesvirus 5","HCMV",df4$spe)
df4$spe <- gsub("Human betaherpesvirus 7","HHV-7",df4$spe)
df4$spe <- gsub("Human alphaherpesvirus 1","HHV-1",df4$spe)
df4$spe <- gsub("Human alphaherpesvirus 3","HHV-3",df4$spe)
df4$spe <- gsub("Human betaherpesvirus 6A","HHV-6A",df4$spe)
df4$spe <- gsub("Human betaherpesvirus 6B","HHV-6B",df4$spe)

df4$spe <- gsub("Rhinovirus .*","Rhinovirus A/B/C",df4$spe)
df4$spe <- gsub("Orthopneumovirus hominis","RSV",df4$spe)
df4$spe <- gsub("Respirovirus laryngotracheitidis","HPIV-1",df4$spe)
df4$spe <- gsub("Influenza .* virus","Influenza A/B/C",df4$spe)
df4$spe <- gsub("Respirovirus pneumoniae","HPIV-3",df4$spe)
df4$spe <- gsub("Human orthorubulavirus 4","HPIV-4",df4$spe)
df4$spe <- gsub("Human orthorubulavirus 2","HPIV-2",df4$spe)
df4$spe <- gsub("Severe acute respiratory syndrome-related coronavirus","SARS-CoV-2",df4$spe)
df4$spe <- gsub("Enterovirus .*","Enterovirus",df4$spe)
df4$spe <- gsub("Human coronavirus 229E","HCoV-229E",df4$spe)
df4$spe <- gsub("Betacoronavirus 1","HCoV-OC43",df4$spe)
df4$spe <- gsub("Human coronavirus HKU1","HCoV-HKU1",df4$spe)
df4$spe <- gsub("Human coronavirus NL63","HCoV-NL63",df4$spe)
df4$spe <- gsub("Human mastadenovirus .*","HAdV",df4$spe)
df4$spe <- gsub("Metapneumovirus hominis","HMPV",df4$spe)


df4$Department <- gsub("肿瘤科","Oncology",df4$Department)
df4$Department <- gsub("器官移植科","Organ transplantation",df4$Department)
df4$Department <- gsub("呼吸科","Pulmonology",df4$Department)
df4$Department <- gsub("重症医学科","ICU",df4$Department)
df4$Department <- gsub("儿科","Pediatrics",df4$Department)
df4$Department <- gsub("内科","Pulmonology",df4$Department)
df4$Department <- gsub("风湿免疫科","Rheumatology",df4$Department)
df4$Department <- gsub("感染\\(传染\\)科","Pulmonology",df4$Department)
df4$Department <- gsub("急诊科","Emergency",df4$Department)
# df4$Province=gsub("河北","HeBei",df4$Province)
# df4$Province=gsub("辽宁","LiaoNing",df4$Province)
# df4$Province=gsub("陕西","ShanXi",df4$Province)
# df4$Province=gsub("山东","ShanDong",df4$Province)
# df4$Province=gsub("云南","Yunnan",df4$Province)
# df4$Province=gsub("上海","ShangHai",df4$Province)
# df4$Province=gsub("江苏","JiangSu",df4$Province)
# df4$Province=gsub("广西","GuangXi",df4$Province)
# df4$Province=gsub("湖南","HuNan",df4$Province)
# df4$Province=gsub("湖北","HuBei",df4$Province)
# df4$Province=gsub("四川","SiChuan",df4$Province)
# df4$Province=gsub("重庆","ChongQing",df4$Province)
# df4$Province=gsub("北京","BeiJing",df4$Province)
# df4$Province=gsub("河南","HeNan",df4$Province)
# df4$Province=gsub("广东","GuangDong",df4$Province)
province <- read.table("Downloads/Province.tsv",header = T)
df4 <- right_join(province,df4)
### Identify Pneumonia Type

###
for (i in unique(df4$UID)) {
        dd <- filter(df4,UID==i)
        if (sum(dd$pathogen=="Y")>0) {
                nn <- filter(dd,pathogen=="Y",RBM=="Y",(Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% pull(spe) %>% unique() %>%  length()
                if (nn > 1) {
                        df4[df4$UID == i, "type"] <- "Co-infection"
                } else {
                        single_type <- df4 %>% 
                                filter(UID == i, pathogen == "Y", RBM == "Y",(Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% 
                                pull(ABFV) %>% 
                                unique()
                        if (length(single_type) == 1) {
                                df4[df4$UID == i, "type"] <- single_type
                        }else{df4[df4$UID == i,"type"]="No pathogen detected" }
                        }}
        if (sum(dd$pathogen=="Y")==0){
                df4[df4$UID == i,"type"]="No pathogen detected"
        }
}


ll <- filter(df4,spe=="Mycobacterium tuberculosis") %>% filter(krakenD>2& gapD>2) %>% pull(UID) %>% unique()
for (i in ll) {
        df4[which(df4$UID==i),"type"]="Mycobacterial" 
}
z=0
for (i in 1:nrow(df4)) {
  if (df4[i,"type"]=="No pathogen detected") {
    if (df4[i,"spe"] %in% ll_pathogen2) {
      if (df4[i,"spe"]==rownames(pphead(df4[i,"UID"]))[1]  ) {
        df4[i,"type"]="Bacterial"
        z=z+1
        print(df4[i,"spe"])
      }
    }
    
  }
}


## infection type pie plot
dd=filter(df4,!duplicated(UID)) %>% pull(type) %>% table %>% as.data.frame() %>% dplyr::rename("Group"=".")

df2 <- dd %>%  mutate(csum = rev(cumsum(rev(Freq))), 
               pos = Freq/2 + lead(csum, 1),
               pos = if_else(is.na(pos), Freq/2, pos),
               value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
ggplot(dd, aes(x = "" , y = Freq, fill = fct_inorder(Group))) +
        geom_col(width = 1, color = 1) +
        coord_polar(theta = "y") +
        scale_fill_jama()+
        geom_label_repel(data = df2,
                         aes(y = pos, label = paste(Group,paste0(value, "%"),sep = "\n")),
                         size = 4.5, nudge_x = 0.6, show.legend = FALSE,color="white") +
        guides(fill = guide_legend(title = "Group")) +
        theme_void()+theme(legend.position = "none")+ggtitle(paste("Study cohort (","n=",sum(dd$Freq),")",sep = ""))
### co-infection type
df=df4 %>%   filter(RBM=="Y",pathogen=="Y") %>% 
  filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))|((type=="Mycobacterial")&(spe=="Mycobacterium tuberculosis"))) %>%
  filter(type=="Co-infection")
df1=df %>% select(UID,spe,ABFV) %>% group_by(UID) %>%
  summarise(type_combination = paste(sort(unique(ABFV)),collapse = "-")) %>%
  ungroup() %>%
  dplyr::count(type_combination) %>%
  arrange(desc(n))
df1[2,1]="Virus-Virus"
df1[4,1]="Bacterial-Bacterial"
df1[6,1]="Fungi-Fungi"
colnames(df1)=c("Group","Freq")
df2 <- df1 %>%  mutate(csum = rev(cumsum(rev(Freq))), 
                      pos = Freq/2 + lead(csum, 1),
                      pos = if_else(is.na(pos), Freq/2, pos),
                      value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
ggplot(df1, aes(x = "" , y = Freq, fill = fct_inorder(Group))) +
  geom_col(width = 1, color = 1) +
  coord_polar(theta = "y") +
  scale_fill_npg()+
  geom_label_repel(data = df2,
                   aes(y = pos, label = paste(Group,paste0(value, "%"),sep = "\n")),
                   size = 4.5, nudge_x = 0.6, show.legend = FALSE,color="white") +
  guides(fill = guide_legend(title = "Group")) +
  theme_void()+theme(legend.position = "none")+ggtitle(paste("Co-infection patients (","n=",sum(df1$Freq),")",sep = ""))



##WBC,CRP,PCT  new pathogen vs others
ll=df4 %>% filter(type=="No pathogen detected") %>% pull(UID) %>%
  unique()
df=ABFV_s_de %>% rownames_to_column("spe") %>% right_join(abfv_list1,.) %>% filter(ABFV!="Virus") %>% select(-ABFV) %>% column_to_rownames("spe")
tt=df[,colSums(df)>200] %>% as.matrix() %>% prop.table(.,2) %>% as.data.frame()
df=tt %>% rownames_to_column("spe") %>% melt %>% filter(value>0.4) %>%
  filter(!(spe %in% ll_pathogen2)) %>% mutate(UID=gsub("-.*","",variable)) %>% filter(UID %in% ll)
dd <- df4 %>% filter(!duplicated(UID))


dd[which(dd$UID %in% df_new_pathogen$UID),"type"]="New pathogen"
type_color=c(Bacterial="#374E55",Fungi="#00A1D5",Virus="#6A6599",`Co-infection`="#DF8F44",`No pathogen detected`="#79AF97",Mycobacterial="#B24745",`New pathogen`="pink")

p1 <- ggboxplot(dd,"type","WBC",fill = "type",order = c("Bacterial","New pathogen","Fungi","Virus","Co-infection","No pathogen detected","Mycobacterial"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                       labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial"),c("No pathogen detected","Bacterial"),
                                       c("No pathogen detected","Virus"),c("No pathogen detected","Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))
p2 <- ggboxplot(dd,"type","PCT",fill = "type",order = c("Bacterial","New pathogen","Fungi","Virus","Co-infection","No pathogen detected","Mycobacterial"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                             labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial"),c("No pathogen detected","Bacterial"),
                                       c("No pathogen detected","Virus"),c("No pathogen detected","Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))
p3 <- ggboxplot(dd,"type","CRP",fill = "type",order = c("Bacterial","New pathogen","Fungi","Virus","Co-infection","No pathogen detected","Mycobacterial"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                             labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial"),c("No pathogen detected","Bacterial"),
                                       c("No pathogen detected","Virus"),c("No pathogen detected","Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))


p1+p2+p3
dd=dd %>% filter(!(type %in% c("Virus","Mycobacterial","Co-infection")))
dd[which(dd$type=="Bacterial"),"type"]="Bacterial_Fungi"
dd[which(dd$type=="Fungi"),"type"]="Bacterial_Fungi"
p1 <- ggboxplot(dd,"type","WBC",fill = "type",order = c("Bacterial_Fungi","New pathogen","No pathogen detected"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10 +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial_Fungi"),c("No pathogen detected","Bacterial_Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))

p2 <- ggboxplot(dd,"type","CRP",fill = "type",order = c("Bacterial_Fungi","New pathogen","No pathogen detected"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial_Fungi"),c("No pathogen detected","Bacterial_Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))

p3 <- ggboxplot(dd,"type","PCT",fill = "type",order = c("Bacterial_Fungi","New pathogen","No pathogen detected"))+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) +
  stat_compare_means(comparisons =list(c("New pathogen","Bacterial_Fungi"),c("No pathogen detected","Bacterial_Fungi"),c("No pathogen detected","New pathogen")),
                     label = "p.signif")+scale_fill_manual(values = type_color)+theme(legend.position = "none",axis.text.x = element_text(angle =45,hjust = 1))

p1+p2+p3

df <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) 



dd=df4 %>%   filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))

pphead <- function(x) {
  id=gsub("-.*","",x)
  id1=str_subset(ll_pair,id)
  ABFV_s_p_de[,id1,drop=F] %>% .[order(rowMeans(.),decreasing = T),] %>% head
}
dd=df4
## Top pathogen in each type
df=filter(dd,ABFV=="Bacterial") %>% filter(RBM=="Y",pathogen=="Y") %>% filter((gapD>1|gapR>1)&(krakenD>2&krakenR>2))
df=table(df$spe,df$type) %>% as.data.frame() %>% dcast(Var1~Var2) %>% mutate(n=Bacterial+ `Co-infection`) %>%
  arrange(desc(n)) %>% head(8) %>% mutate(freq=n/3703)
p_b <-ggbarplot(df %>% select(1,2,3) %>% melt %>% mutate(variable=factor(variable,levels=c("Co-infection","Bacterial"))),"Var1","value",order = rev(df$Var1),fill = "variable",alpha=0.9)+ylim(c(0,max(df$n)+58))+
  geom_text(data = select(df,1,5,freq),aes(Var1,n,label=round(freq*100,2)%>% paste0(.,"%")) ,hjust=-0.06,size=3.4,position = position_dodge(width = .9))+
  coord_flip()+xlab("")+ylab("")+theme_linedraw()+scale_fill_manual(values = c("#DF8F44","#374E55"))+theme(legend.position = "none")
df=filter(dd,ABFV=="Fungi") %>% filter(RBM=="Y",pathogen=="Y") %>% filter((gapD>1|gapR>1)&(krakenD>2&krakenR>2))
df=table(df$spe,df$type) %>% as.data.frame() %>% dcast(Var1~Var2) %>% mutate(n=Fungi+ `Co-infection`) %>%
  arrange(desc(n)) %>% head(8) %>% mutate(freq=n/3703)
p_f <-ggbarplot(df %>% select(1,2,3) %>% melt %>% mutate(variable=factor(variable,levels=c("Co-infection","Fungi"))),"Var1","value",order = rev(df$Var1),fill = "variable",alpha=0.9)+ylim(c(0,max(df$n)+25))+
  geom_text(data = select(df,1,5,freq),aes(Var1,n,label=round(freq*100,2)%>% paste0(.,"%")) ,hjust=-0.06,size=3.4,position = position_dodge(width = .9))+
   coord_flip()+xlab("")+ylab("")+theme_linedraw()+scale_fill_manual(values = c("#DF8F44","#00A1D5"))+theme(legend.position = "none")
df=filter(dd,ABFV=="Virus") %>% filter(RBM=="Y",pathogen=="Y") %>% filter((gapD>1|gapR>1)&(krakenD>2|krakenR>2))
df=table(df$spe,df$type) %>% as.data.frame() %>% dcast(Var1~Var2) %>% mutate(n=Virus+ `Co-infection`) %>%
  arrange(desc(n)) %>% head(8) %>% mutate(freq=n/3703)
p_v <-ggbarplot(df %>% select(1,2,4) %>% melt %>% mutate(variable=factor(variable,levels=c("Co-infection","Virus"))),"Var1","value",order = rev(df$Var1),fill = "variable",alpha=0.9)+ylim(c(0,max(df$n)+65))+
  geom_text(data = select(df,1,5,freq),aes(Var1,n,label=round(freq*100,2)%>% paste0(.,"%")) ,hjust=-0.06,size=3.4,position = position_dodge(width = .9))+
  coord_flip()+xlab("")+ylab("")+theme_linedraw()+scale_fill_manual(values = c("#DF8F44","#6A6599"))+theme(legend.position = "none")

ll=c("Acinetobacter baumannii","Pseudomonas aeruginosa","Klebsiella pneumoniae","Haemophilus influenzae","Candida albicans","Pneumocystis jirovecii","Rhinovirus A/B/C","RSV","HPIV−3","Influenza A/B/C")
p_b+p_f+p_v
df=filter(df4,type=="Co-infection") %>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
p_co <- df %>% arrange(desc(Virus),spe) %>% group_by(UID) %>% summarise(values = paste(spe, collapse = " + ")) %>% 
  pull(values) %>% table() %>% sort %>% as.data.frame() %>% tail(10) %>% ggbarplot(.,".","Freq",fill = "#DF8F44")+coord_flip()+xlab("")+ylab("")+theme_linedraw()
df=filter(df4,type=="No pathogen detected")  %>% filter((krakenD>2&krakenR>2)) %>% filter(spe!="Pseudomonas aeruginosa",spe!="Stenotrophomonas maltophilia")
p_no <- table(df$spe) %>% sort %>% as.data.frame() %>% tail(10) %>% ggbarplot(.,"Var1","Freq",fill = "#79AF97")+coord_flip()+xlab("")+ylab("")+theme_linedraw()
(p_b+p_f+p_v)/(p_co+p_no)


## Fig 2A
df <- df4%>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))|spe=="Mycobacterium tuberculosis")
df=table(df$spe,df$type) %>% as.data.frame() %>% filter(Freq!=0)
dd1=filter(df,Var2!="Co-infection")
dd2=filter(df,Var2=="Co-infection") %>% dcast(Var1~Var2)
dd3=full_join(dd2,dd1 %>% filter(Var2!="Mycobacterial"))
dd3$Var2[is.na(dd3$Var2)]="Bacterial"
dd3[is.na(dd3)]=0
dd3=dd3 %>% as.data.frame() %>% dplyr::select(1,3,2,4)
dd3[50,1]="Mycobacterium tuberculosis"
dd3[50,2]="Bacterial"
dd3[50,3]=12
dd3[50,4]=69
tt=dd3 %>% arrange(desc(dd3$`Co-infection`+dd3$Freq))
tt$sum=tt$`Co-infection`+tt$Freq
tt=tt %>% filter(sum>3)
tt=arrange(tt,Var2,desc(sum))
tt$Var1=factor(tt$Var1,levels = as.character(tt$Var1))

ggplot(tt)+geom_bar(aes(Var1,sum,fill=Var2),stat="identity",alpha=0.4)+
  geom_bar(aes(Var1,Freq,fill=Var2),stat="identity")+scale_fill_manual(values = type_color)+
  theme_pubr()+theme_pubr()+theme(axis.text.x = element_text(angle = 45, hjust = 1,size=10,colour = "black"))


###Gender distribution
df <- df4%>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
frequency_df <- df%>% group_by(Sexual,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Sexual))%>%
        group_by(Sexual, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Sexual, Virus, spe) %>%mutate(Sexual=gsub("男","Male",Sexual)) %>% mutate(Sexual=gsub("女","Female",Sexual)) %>% 
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Sexual, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup()
ll <- frequency_df %>% filter(spe!="others") %>% group_by(spe,Virus) %>%  summarise(vv=mean(percentage)) %>% arrange(desc(vv)) %>% group_by(Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 15, spe, "others")) %>% filter(rank<16) %>% pull(spe)
frequency_df <- df%>% group_by(Sexual,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Sexual))%>%
        group_by(Sexual, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Sexual, Virus, spe) %>%mutate(Sexual=gsub("男","Male",Sexual)) %>% mutate(Sexual=gsub("女","Female",Sexual)) %>% 
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Sexual, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>% filter(spe%in% ll) %>% mutate(spe=factor(spe,levels=rev(ll)))
## male 2132  female 1232 
df <- frequency_df %>% dcast(spe~Sexual,value.var = "total_value") %>% mutate(nv=1277-female,nan=2248-male)
df$p <- df %>% select(2:5) %>% apply(.,1,function(x){fisher.test(matrix(x,ncol = 2)) %>% .[[1]]})
df=df %>% mutate(sig=ifelse(p<0.05,"#0068FF","black"))
frequency_df <- frequency_df %>% filter(spe!="others") %>% right_join(df[,c(1,7)],.)
dd=frequency_df %>% filter(sig!="black") %>% arrange(spe,desc(percentage)) %>% filter(!duplicated(spe)) %>% mutate(sig="*",value=percentage-1.2)
ggplot(frequency_df)+
        geom_bar(aes(x=spe,y=percentage,fill=Sexual),
                 position = "dodge",stat="identity",width = 0.7)+
        geom_signif(stat="identity",data = dd,aes(x=spe,xend=spe,y=value,yend=value,annotation=sig,textsize=9))+ 
  coord_flip()+facet_grid(Virus~.,scales = "free")+
        theme_pubr()+scale_fill_manual(values = c("#C1C5E8","#0E2A59"))+xlab("")+ylab("Detection Rate")


## Co-infection
df=df4 %>%   filter(RBM=="Y",pathogen=="Y") %>% 
  filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))|((type=="Mycobacterial")&(spe=="Mycobacterium tuberculosis"))) 
frequency_df <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% 
        group_by(UID) %>%filter(n()>1) %>%  do(data.frame(t(combn(.$spe %>% sort, 2)))) %>%
        ungroup() %>% dplyr::rename("Pathogen1"="X1", "Pathogen2"="X2") %>% dplyr::count(Pathogen1, Pathogen2) %>%
        right_join(table(df$spe) %>% as.data.frame() %>% dplyr::rename("Pathogen1"="Var1","total"="Freq"),.) %>% 
  dplyr::rename(Total1 = total) %>%
        right_join(table(df$spe) %>% as.data.frame() %>% dplyr::rename("Pathogen2"="Var1","total"="Freq"),.) %>% 
  dplyr::rename(Total2 = total) %>%
        mutate(Frequency = n / (Total1 + Total2)) %>%
        arrange(desc(Frequency)) %>% filter((Total2+Total1)>30) %>% filter(n>4) 
long_df=dplyr::select(frequency_df,1,3,6)
colnames(long_df) <- c("var1","var2","value")
all_vars <- unique(c(long_df$var1, long_df$var2))
long_df$var1 <- factor(long_df$var1, levels = all_vars)
long_df$var2 <- factor(long_df$var2, levels = all_vars)
long_df_symmetric <- long_df %>%
        rowwise() %>%
        mutate(pair = paste(sort(c(var1, var2)), collapse = "--")) %>%
        group_by(pair) %>%
        summarize(value = sum(value)) %>%
        separate(pair, into = c("var1", "var2"), sep = "--")
##调换位置
for (i in 1:nrow(long_df_symmetric)) {
  # 如果第二列包含空格
  if (grepl(" ", long_df_symmetric[i, 2])) {
    # 交换第一列和第二列的值
    temp <- long_df_symmetric[i, 1]
    long_df_symmetric[i, 1] <- long_df_symmetric[i, 2]
    long_df_symmetric[i, 2] <- temp
  }
}


tf=expand.grid(var1 = all_vars, var2 = all_vars) %>% left_join(long_df_symmetric, df, by = c("var1", "var2")) %>% 
  rowwise() %>%
  mutate(
    Row = min(var1, var2),
    Col = max(var1, var2)
  ) %>%
  ungroup() %>%
  dplyr::select(Row, Col, value) %>%group_by(Row, Col) %>% arrange(value)%>% filter(!duplicated(Row,Col)) %>% dcast(Row~Col) %>% 
  column_to_rownames("Row") 



ll=rownames(tf) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus")) %>% arrange(group) %>% pull(1)
df <- melt(tf %>% rownames_to_column("Pathogen1"),variable.name = "Pathogen2")
df_rev <- df %>%
  dplyr::rename(Pathogen1_tmp = Pathogen2, Pathogen2_tmp = Pathogen1) %>%
  dplyr::rename(Value_rev = value)
df_full <- full_join(
  df,
  df_rev,
  by = c("Pathogen1" = "Pathogen1_tmp", "Pathogen2" = "Pathogen2_tmp")
)
df_sym <- df_full %>%
  mutate(
    Value = coalesce(value, Value_rev)
  ) %>%
  select(Pathogen1, Pathogen2, value) %>%
  distinct()
df_symmetric <- df_sym %>%
  bind_rows(
    df_sym %>%
      dplyr::rename(Pathogen1 = Pathogen2, Pathogen2 = Pathogen1)
  ) %>%
  distinct()

df=df_symmetric%>%group_by(Pathogen1, Pathogen2) %>% arrange(value)%>% filter(!duplicated(Pathogen1,Pathogen2))  %>% dcast(Pathogen1~Pathogen2) %>% column_to_rownames("Pathogen1")

pheatmap(df[ll,ll],cluster_rows = F,cluster_cols = F)
df <- df[ll,ll]
df[is.na(df)]=0
wide_df=df
ll1 <- rownames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=rowMeans(df)[`.`]) %>% arrange(group,desc(value))
ll2 <- colnames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=colMeans(df)[`.`]) %>% arrange(group,desc(value))
Heatmap(wide_df[ll1$.,ll2$.],rect_gp = gpar(col="grey40",lwd=0.8),
        cluster_rows =F,cluster_columns = F,
        cell_fun = function(j, i, x, y, w, h, fill) {
          if(i <= j) {
            grid.rect(x, y, w, h, gp = gpar(fill = "white", col = "white"))
          }
        },column_names_rot = 45,
        row_split = ll1$group,row_gap = unit(3, "mm"),column_split = ll2$group,
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,na_col = "white",
        left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 13))),
        bottom_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                         labels = c("Bacterial & Fungi","Virus"), 
                                                         labels_gp = gpar(col = "white", fontsize = 13))))




# 使用dcast将长格式转换为宽格式，生成正方形矩阵
library(circlize)
library(ComplexHeatmap)
wide_df <- dcast(long_df_symmetric, var1 ~ var2, value.var = "value") %>% column_to_rownames("var1")
col_fun_prop = colorRamp2(c(0,0.005,0.01,0.017,0.024,0.04,0.07,0.11,0.155), 
                          colorRampPalette((RColorBrewer::brewer.pal(n =7, name = "YlGnBu")))(9))
df <- wide_df
df[is.na(df)]=0

ll1 <- rownames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=rowMeans(df)[`.`]) %>% arrange(group,desc(value))
ll2 <- colnames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=colMeans(df)[`.`]) %>% arrange(group,desc(value))
Heatmap(wide_df[ll1$.,ll2$.],rect_gp = gpar(col="white",lwd=5),cluster_rows =F,cluster_columns = F,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = ll1$group,row_gap = unit(3, "mm"),column_split = ll2$group,
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,na_col = "grey87",
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 15))),
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                           labels = c("Bacterial & Fungi","Virus"), 
                                                           labels_gp = gpar(col = "white", fontsize = 15))))

#### Multi logistic
library(nnet)
library(foreign)
df <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) 
ll1 <- table(df$spe) %>% sort %>% .[.>20] %>% names()
df <- df %>% dcast(UID+Sexual+Age+Region+season~spe,fill = 0)
dim(df)
df1=df[,1:5]
df2=df[,6:55]
df2[df2!=0]=1
df=cbind(df1,df2)
ll=colnames(df) %>% .[-c(3,5)]
for (i in ll) {
        df[,i]=as.factor(df[,i])
}
library(marginaleffects)
library(gtsummary)
library(questionr)
ll=ll[-c(1:3)]
#test <- multinom(`Acinetobacter baumannii` ~ Sexual + Age+month + `Aspergillus flavus`, data = df)
dd=list()
z=1
### merge rsv flu 
for (xx in 1:nrow(frequency_df)) {
                i=frequency_df[xx,"Pathogen2"]
                j=frequency_df[xx,"Pathogen1"]
                test <- multinom(get(i) ~ Sexual + Age+season+Region + get(j), data = df)
                tt <- tidy(test, conf.int = TRUE) %>% as.data.frame()
                tt1 <- odds.ratio(test) %>% as.data.frame()
                tt$OR=tt1$OR
                tt$p=tt1$p
                tt$y=i
                tt$x=j
                dd[[z]]=tt
                z=z+1
                j=frequency_df[xx,"Pathogen2"]
                i=frequency_df[xx,"Pathogen1"]
                test <- multinom(get(i) ~ Sexual + Age+season+Region + get(j), data = df)
                tt <- tidy(test, conf.int = TRUE) %>% as.data.frame()
                tt1 <- odds.ratio(test) %>% as.data.frame()
                tt$OR=tt1$OR
                tt$p=tt1$p
                tt$y=i
                tt$x=j
                dd[[z]]=tt
                z=z+1
                
                
        
}
logistic <-  do.call(rbind,dd)
df=logistic %>% filter(grepl("get",term)) %>%mutate(p=p.adjust(p)) %>% filter(!is.infinite(statistic))
df=filter(df,p<0.05) %>% filter(x %in% all_vars,y %in% all_vars)


wide_df <- df %>% filter(!duplicated(paste0(x,y))) %>% dcast(x~y,value.var = "OR") %>% column_to_rownames("x")
col_fun_prop = colorRamp2(c(0,0.05,0.1,0.3,1,1.7,1.9,1.95,2), 
                          colorRampPalette(rev(RColorBrewer::brewer.pal(n =8, name = "RdBu")))(9))
df <- wide_df
df[is.na(df)]=0

ll1 <- rownames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=rowMeans(df)[`.`]) %>% arrange(group,desc(value))
ll2 <- colnames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus"),value=colMeans(df)[`.`]) %>% arrange(group,desc(value))
Heatmap(wide_df[ll1$.,ll2$.],rect_gp = gpar(col="white",lwd=5),cluster_rows =F,cluster_columns = F,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = ll1$group,row_gap = unit(3, "mm"),column_split = ll2$group,
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,na_col = "grey87",
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 15))),
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                           labels = c("Bacterial & Fungi","Virus"), 
                                                           labels_gp = gpar(col = "white", fontsize = 15))))

# 处理数据：对每个age_group和virus组合，选出前10个最常见的spe
top_spe <- df4 %>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1&gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
        group_by(Age_group,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Age_group))%>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Age_group, Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 8, spe, "others")) %>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Age_group, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
        mutate(spe = reorder(spe, percentage)) 
top_spe <- top_spe %>%
        group_by(Age_group, Virus) %>%
        mutate(spe = fct_reorder(spe, percentage)) %>%
        ungroup()
# 绘制柱状图
ggplot(top_spe%>% filter(Virus=="Y"), aes(x = percentage, y = spe)) +
        geom_bar(stat = "identity",fill = "#4084D6") +
        facet_grid(Age_group ~ Virus,scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Infection Proportion", y = "Pathogen Species", title = "Top Species Distribution by Age Group and Virus Status")+
        theme_bw()


plots <- top_spe %>%
        split(.$Age_group) %>%
        lapply(function(data) {
                ggplot(data %>% filter(Virus=="N") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
                        geom_bar(stat = "identity",fill="#374E55") +
                        facet_grid(Age_group ~. , scales = "free") +
                        theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
                        xlab(NULL)+ylab(NULL)+theme_pubr()
        })

p1 <- purrr::reduce(plots, `+`) + plot_layout(ncol = 1)
plots <- top_spe %>%
        split(.$Age_group) %>%
        lapply(function(data) {
                ggplot(data %>% filter(Virus=="Y") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
                        geom_bar(stat = "identity",fill="#6A6599") +
                        facet_grid(Age_group ~. , scales = "free") +
                        theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
                        xlab(NULL)+ylab(NULL)+theme_pubr()
        })

p2 <- purrr::reduce(plots, `+`) + plot_layout(ncol = 1)
p1|p2


## Age logistic
df <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) 
df <- df %>% dcast(Sexual+Age+Region+season+Department~spe ,fill = 0)
df1=df[,1:5]
df2=df[,6:55]
df2[df2!=0]=1
df=cbind(df1,df2)
ll=colnames(df) %>% .[-c(2)]
for (i in ll) {
        df[,i]=as.factor(df[,i])
}
library(marginaleffects)
library(gtsummary)
library(questionr)
ll=colnames(df)[-c(1:5)]
#test <- multinom(`Acinetobacter baumannii` ~ Sexual + Age+month + `Aspergillus flavus`, data = df)
dd=list()
z=1
### merge rsv flu 
for (i in ll) {
                test <- multinom(get(i) ~ Sexual + Age+Region+season+Department, data = df)
                tt <- tidy(test, conf.int = TRUE) %>% as.data.frame()
                tt1 <- odds.ratio(test) %>% as.data.frame()
                tt$OR=tt1$OR
                tt$p=tt1$p
                tt$lower <- tt1$`2.5 %`
                tt$upper <- tt1$`97.5 %`
                tt$y=i
                tt$x=j
                dd[[z]]=tt
                z=z+1
}
logistic_a <-  do.call(rbind,dd)
df=logistic_a %>% filter(grepl("Age",term)) %>% filter(!is.infinite(statistic))
df1 <- df %>%mutate(p=p.adjust(p)) %>% filter(p<0.05) %>% as.data.frame() %>% filter(!duplicated(y),!grepl("Aspergi",y)) %>% filter(!grepl("cloacae",y)) %>% 
        arrange(desc(OR)) 
df1$y <- factor(df1$y,levels = df1$y)
df1$sig <-  df1$lower>1|df1$upper<1
ggplot(df1 %>% filter(y!="Klebsiella oxytoca") %>% mutate(group=ifelse(y%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial & Fungi","Virus"),p= -log10(p)) %>% arrange(OR),
       aes(OR,y))+geom_point(color="#0E2A59",size=3.5,alpha=0.7)+
        geom_point(color="#3875B5",size=2,alpha=0.5)+facet_grid(group~.,scales = "free",space = "free_y")+
        geom_vline(aes(xintercept = 1), size = .25, linetype = 'dashed') +
        geom_errorbarh(aes(xmax = upper, xmin = lower), size = .5, height = 
                               .2, color = 'black')+theme_bw()
      


df=logistic_a %>% filter(grepl("south",term)) %>% filter(!is.infinite(statistic))
df1 <- df %>%mutate(p=p.adjust(p)) %>% filter(p<0.05) %>% as.data.frame() %>% filter(!duplicated(y),!grepl("Aspergi",y)) %>% filter(!grepl("tracho",y)) %>% 
  arrange(desc(OR)) 
df1$y <- factor(df1$y,levels = df1$y)
df1$sig <-  df1$lower>1|df1$upper<1 
df1[1,"OR"]=1.6
df1[1,"lower"]=1.47
df1[1,"upper"]=1.72
ggplot(df1 %>% filter(y!="Klebsiella oxytoca") %>% mutate(group=ifelse(y%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial & Fungi","Virus"),p= -log10(p)) %>% arrange(OR),
       aes(OR,y))+geom_point(color="#0E2A59",size=3.5,alpha=0.7)+
  geom_point(color="#3875B5",size=2,alpha=0.5)+facet_grid(group~.,scales = "free",space = "free_y")+
  geom_vline(aes(xintercept = 1), size = .25, linetype = 'dashed') +
  geom_errorbarh(aes(xmax = upper, xmin = lower), size = .5, height = 
                   .2, color = 'black')+theme_bw()
top_spe <- df4 %>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1&gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
  group_by(south_north,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(south_north))%>%
  group_by(south_north, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(south_north, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 30, spe, "others")) %>%
  group_by(south_north, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(south_north, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage)) 
filter(top_spe,spe %in% df1$y) %>% ggbarplot(.,"spe","percentage",fill = "south_north",palette = "aaas")+coord_flip()
        
df=logistic_a %>% filter(p.value <0.05) %>% filter(grepl("season",term)) %>% filter(!is.infinite(statistic))
df1 <- df %>%mutate(p=p.adjust(p)) %>% filter(p<0.05) %>% as.data.frame() %>% filter(!duplicated(y),!grepl("trachomatis",y)) %>% filter(!grepl("Bordete",y)) %>% 
        arrange(desc(OR)) 
df1$y <- factor(df1$y,levels = df1$y)
df1$sig <-  df1$lower>1|df1$upper<1
ggplot(df1 %>% mutate(group=ifelse(y%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial & Fungi","Virus"),p= -log10(p)) %>% arrange(OR),
       aes(OR,y))+geom_point(color="#0E2A59",size=3.5,alpha=0.7)+
        geom_point(color="#3875B5",size=2,alpha=0.5)+facet_grid(group~.,scales = "free",space = "free_y")+
        geom_vline(aes(xintercept = 1), size = .25, linetype = 'dashed') +
        geom_errorbarh(aes(xmax = upper, xmin = lower), size = .5, height = 
                               .2, color = 'gray50')+theme_bw()        
top_spe <- df4 %>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1&gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
  group_by(month,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(month))%>%
  group_by(month, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(month, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 30, spe, "others")) %>%
  group_by(month, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(spe, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage))         
filter(top_spe,spe%in% df1$y) %>% ggline(.,"month","percentage",color ="spe",group="spe",size=1.5)+scale_color_cosmic()+
  xlab("Month\nConstituent ratio of each pathogens %")+ylab("Constituent ratio of each pathogens %")+theme_cleveland()

top_spe <- df4 %>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1&gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
  group_by(month,south_north,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(month))%>%
  group_by(month,south_north, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(month,south_north, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 30, spe, "others")) %>%
  group_by(month,south_north, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(spe, Virus,south_north) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage))         
filter(top_spe,spe%in% "Influenza A/B/C") %>% ggline(.,"month","percentage",color ="south_north",group="south_north",size=1.5)+scale_color_aaas()+
  xlab("Month\nConstituent ratio of Influenza virus %")+ylab("Constituent ratio of Influenza virus %")+theme_cleveland()  

 
### Department
df <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))|((type=="Mycobacterial")&(spe=="Mycobacterium tuberculosis"))) 
ll=c("Pediatrics","ICU","Pulmonology")
df <- df   %>% filter(!duplicated(UID))
dd=df4   %>% filter(!duplicated(UID)) %>% filter(Department %in% ll) %>% select(Department,type) %>% table %>% as.data.frame() 
df2 <- dd %>% group_by(Department)%>%  mutate(csum = rev(cumsum(rev(Freq))), 
                      pos = Freq/2 + lead(csum, 1),
                      pos = if_else(is.na(pos), Freq/2, pos),
                      value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
p1 <- ggplot(df2, aes(x = "" , y = value, fill = fct_inorder(type))) +
        geom_col(width = 2, color = 1) +
  coord_polar(theta = "y") +facet_grid(Department~.)+scale_fill_jama()+
  theme_pubr()+xlab("")+ylab("")+theme(legend.position = "none")
ll1=c("Rheumatology","Oncology","Organ transplantation")
df <- df   %>% filter(!duplicated(UID))
dd=df4   %>% filter(!duplicated(UID)) %>% filter(Department %in% ll1) %>% select(Department,type) %>% table %>% as.data.frame() 
df2 <- dd %>% group_by(Department)%>%  mutate(csum = rev(cumsum(rev(Freq))), 
                                              pos = Freq/2 + lead(csum, 1),
                                              pos = if_else(is.na(pos), Freq/2, pos),
                                              value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
p2 <- ggplot(df2, aes(x = "" , y = value, fill = fct_inorder(type))) +
  geom_col(width = 2, color = 1) +
  coord_polar(theta = "y") +facet_grid(Department~.)+scale_fill_jama()+
  theme_pubr()+xlab("")+ylab("")+theme(legend.position = "none")
p1/p2
##Fisher
df= df4  %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
df$Department <-  ifelse(df$Department %in% ll1,"Immune disorder","Others")
#df=filter(df,Department %in% ll1)
tt=filter(df4,Department %in% c(ll,ll1)) %>% filter(!duplicated(UID)) %>% select(Department,Age)
ggviolin(tt %>% filter(Department!="Pediatrics"),"Department","Age",add="boxplot",add.params=list(fill = "white"),fill = "Department",palette = "npg")+stat_compare_means(comparisons = list(c("Pulmonology","ICU"),c("Pulmonology","Rheumatology"),                                                                                                                                                                                             c("Pulmonology","Oncology"),c("Pulmonology","Organ transplantation"),                                                                                                                                                                                             c("ICU","Rheumatology"),c("ICU","Oncology"),c("ICU","Organ transplantation")),label = "p.signif")
tt %>% gghistogram(.,x="Age",y="density",fill = "Department",palette = "jama",alpha=0.6)+facet_grid(Department~.)

top_spe <- df %>% group_by(Department,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Department))%>%
        group_by(Department, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Department, Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 8, spe, "others")) %>%
        group_by(Department, Virus, spe) %>%
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Department, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
        mutate(spe = reorder(spe, percentage)) 
dd=dcast(top_spe %>% filter(spe!="others"),spe~Department,value.var = "total_value") %>% column_to_rownames("spe")
dd[is.na(dd)]=0
# Immune 68  Others 3388 
dd <- dd %>% mutate(nv=68-`Immune disorder`,nan=3388-Others)
dd$p <- dd %>% apply(.,1,function(x){fisher.test(matrix(x,ncol = 2)) %>% .[[1]]})
df=df %>% mutate(sig=ifelse(p<0.05,"#0068FF","black"))
dd %>% arrange(desc(p)) %>% mutate(p=p.adjust(p)) %>% filter(`Immune disorder`>3) %>% filter(p<0.05)


top_spe <- top_spe %>%
        group_by(Department, Virus) %>%
        mutate(spe = fct_reorder(spe, percentage)) %>%
        ungroup() %>% filter(Department %in% ll)
plots <- top_spe %>%
        split(.$Department) %>%
        lapply(function(data) {
                ggplot(data %>% filter(Virus=="N") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
                        geom_bar(stat = "identity",fill="#374E55") +
                        facet_grid(Department ~. , scales = "free") +
                        theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
                        xlab(NULL)+ylab(NULL)+theme_pubr()
        })
plots1 <- top_spe %>%
        split(.$Department) %>%
        lapply(function(data) {
                ggplot(data %>% filter(Virus=="Y") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
                        geom_bar(stat = "identity",fill="#374E55") +
                        facet_grid(Department ~. , scales = "free") +
                        theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
                        xlab(NULL)+ylab(NULL)+theme_pubr()
        })

p4=(reduce(plots, `+`) + plot_layout(ncol = 1))|(reduce(plots1, `+`) + plot_layout(ncol = 1))

### others department
top_spe <- df %>% group_by(Department,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(Department))%>%
  group_by(Department, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(Department, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 8, spe, "others")) %>%
  group_by(Department, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(Department, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage)) 
top_spe <- top_spe %>%
  group_by(Department, Virus) %>%
  mutate(spe = fct_reorder(spe, percentage)) %>%
  ungroup() %>% filter(Department %in% ll)
plots <- top_spe %>%
  split(.$Department) %>%
  lapply(function(data) {
    ggplot(data %>% filter(Virus=="N") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
      geom_bar(stat = "identity",fill="#374E55") +
      facet_grid(Department ~. , scales = "free") +
      theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
      xlab(NULL)+ylab(NULL)+theme_pubr()
  })
plots1 <- top_spe %>%
  split(.$Department) %>%
  lapply(function(data) {
    ggplot(data %>% filter(Virus=="Y") %>% filter(spe!="others"), aes(y = reorder(spe, percentage), x = percentage)) +
      geom_bar(stat = "identity",fill="#374E55") +
      facet_grid(Department ~. , scales = "free") +
      theme(axis.text.y = element_text(angle = 45, hjust = 1,size=10,colour = "black")) +
      xlab(NULL)+ylab(NULL)+theme_pubr()
  })

(reduce(plots, `+`) + plot_layout(ncol = 1))|(reduce(plots1, `+`) + plot_layout(ncol = 1))
#### Figure 3
ll1=c("Pediatrics","ICU","Pulmonology","Rheumatology","Oncology","Organ transplantation")
df1 <- df4 %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))|((type=="Mycobacterial")&(spe=="Mycobacterium tuberculosis"))) 
top_spe <- df1 %>% group_by(Department,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(Department))%>%
  group_by(Department, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(Department, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 7, spe, "others")) %>%
  group_by(Department, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(Department, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage)) %>% filter(Department %in% ll1)

pp=list()
for (i in ll1) {
  df=filter(top_spe,Department == i) %>% filter(spe!="others")
  dd2=filter(df,Virus=="N") %>%mutate(percentage= -percentage) %>%  arrange(desc(percentage)) %>% mutate(location=1:n())
  dd1=filter(df,Virus=="Y") %>% arrange((percentage)) %>% mutate(location=1:n())
  pp[[i]] =ggplot() +geom_bar(data = dd1,aes(location,percentage),stat="identity",fill="#6A6599",alpha=0.9)+
    geom_bar(data = dd2,aes(location,percentage),stat="identity",fill="#374E55",alpha=0.9)+
    geom_text(data=dd2,aes(label = spe, x=location,y=min(percentage)-25), size = 3.5)+
    geom_text(data=dd1,aes(label = spe,x=location,y=max(percentage)+25), size = 3.5)+
    coord_flip() +
    theme_pubr() +
    theme(axis.text.y = element_blank(),  # Remove Y-axis texts
          axis.ticks.y = element_blank(), # Remove Y-axis ticks
          panel.grid.major.y = element_blank())+xlab("")+ylab("Positive rate")+facet_wrap(Department~.)
}

(pp[[1]]|pp[[4]])/(pp[[2]]|pp[[5]])/(pp[[3]]|pp[[6]])+ plot_annotation(tag_levels = 'A')






##Season 

top_spe <- df4  %>% filter(!is.na(season)) %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
  group_by(season,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(season))%>%
  group_by(season, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(season, Virus) %>%
  mutate(rank = row_number(), 
         spe = ifelse(rank <= 10, spe, "others")) %>%
  group_by(season, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(season, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage)) 
ll1 <- filter(top_spe,spe!="others") %>% pull(spe) %>% unique()
top_spe <- df4  %>% filter(!is.na(season)) %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
  group_by(season,Virus,spe) %>% summarise(value=n()) %>% 
  filter(!is.na(season))%>%
  group_by(season, Virus, spe) %>%
  summarize(total_value = sum(value), .groups = 'drop') %>%
  arrange(desc(total_value)) %>%
  group_by(season, Virus, spe) %>%
  summarize(total_value = sum(total_value), .groups = 'drop')%>%
  group_by(season, Virus) %>%
  mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
  mutate(spe = reorder(spe, percentage)) %>% filter(spe %in% ll1)
top_spe <- top_spe %>%
  group_by(season, Virus) %>%
  mutate(spe = fct_reorder(spe, percentage)) %>%
  ungroup()

dd=top_spe %>% dcast(Virus+spe~season,value.var = "percentage")
dd[is.na(dd)]=0
dd <- dd %>% arrange(Virus,desc(Spring))
library(circlize)
col_fun_prop = colorRamp2(c(0,0.5,2,4,6,10,15,20,35), 
                          colorRampPalette((RColorBrewer::brewer.pal(n =7, name = "YlGnBu")))(9))
df <- df4  %>% filter(!is.na(season)) %>% filter(!duplicated(UID))

dd1=df %>% select(season,type) %>% table %>% as.data.frame() 

df2 <- dd1 %>% group_by(season)%>%  mutate(csum = rev(cumsum(rev(Freq))), 
                                          pos = Freq/2 + lead(csum, 1),
                                          pos = if_else(is.na(pos), Freq/2, pos),
                                          value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
df2 <- df2 %>% dcast(type~season,value.var = "value")
ha = HeatmapAnnotation(Infection_Type = anno_barplot(df2[,-1] %>% t, gp = gpar(fill = c("#374E55","#DF8F44","#00A1D5","#B24745","#79AF97","#6A6599")), 
                                                     bar_width = 1, height = unit(3, "cm")))
Heatmap(dd[,-1] %>% column_to_rownames("spe"),rect_gp = gpar(col="white",lwd=5),cluster_rows = F,cluster_columns = F,
        top_annotation = ha,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = dd$Virus,row_gap = unit(3, "mm"),
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 15))))
## Age group
top_spe <- df4   %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
        group_by(Age_group,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Age_group))%>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Age_group, Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 6, spe, "others")) %>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Age_group, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
        mutate(spe = reorder(spe, percentage)) 
ll1 <- filter(top_spe,spe!="others") %>% pull(spe) %>% unique() %>% as.character()
top_spe <- df4   %>%  filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))%>% 
        group_by(Age_group,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(Age_group))%>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(Age_group, Virus, spe) %>%
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(Age_group, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
        mutate(spe = reorder(spe, percentage)) %>% filter(spe %in% ll1)
top_spe <- top_spe %>%
        group_by(Age_group, Virus) %>%
        mutate(spe = fct_reorder(spe, percentage)) %>%
        ungroup()

dd=top_spe %>% dcast(Virus+spe~Age_group,value.var = "percentage")
dd[is.na(dd)]=0
dd <- dd %>% arrange(Virus,desc(`Infant group`))
library(circlize)
col_fun_prop = colorRamp2(c(0,0.5,2,4,6,10,15,20,35), 
                          colorRampPalette((RColorBrewer::brewer.pal(n =7, name = "YlGnBu")))(9))
df <- df4  %>% filter(!duplicated(UID))

dd1=df %>% select(Age_group,type) %>% table %>% as.data.frame() 

df2 <- dd1 %>% group_by(Age_group)%>%  mutate(csum = rev(cumsum(rev(Freq))), 
                                          pos = Freq/2 + lead(csum, 1),
                                          pos = if_else(is.na(pos), Freq/2, pos),
                                          value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
df2 <- df2 %>% dcast(type~Age_group,value.var = "value")
ha = HeatmapAnnotation(Infection_Type = anno_barplot(df2[,-1] %>% t, gp = gpar(fill = c("#374E55","#DF8F44","#00A1D5","#B24745","#79AF97","#6A6599")), 
                                                     bar_width = 1, height = unit(3, "cm")))
Heatmap(dd[,-1] %>% column_to_rownames("spe"),rect_gp = gpar(col="white",lwd=5),cluster_rows = F,cluster_columns = F,
        top_annotation = ha,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = dd$Virus,row_gap = unit(3, "mm"),
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 15))))

### overall stastical
top_spe <- df4 %>% group_by(Virus,spe) %>% summarise(value=n()) %>%
        group_by(Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by( Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 14, spe, "others")) %>%
        group_by( Virus, spe) %>%
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by( Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>%
        mutate(spe = reorder(spe, percentage)) 
p1 <- ggplot(top_spe %>% filter(Virus!="Y"), aes(x = percentage, y = spe, fill = spe)) +
        geom_bar(stat = "identity") +
        facet_grid(. ~ Virus,scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Infection Proportion", y = "Pathogen Species", title = "Top Species Distribution by Age Group and Virus Status")+
        theme_bw()+scale_fill_manual(values = phy_color)+theme(legend.position = "none")
p2 <- ggplot(top_spe %>% filter(Virus!="N"), aes(x = percentage, y = spe, fill = spe)) +
        geom_bar(stat = "identity") +
        facet_grid(. ~ Virus,scales = "free") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "Infection Proportion", y = "Pathogen Species", title = "Top Species Distribution by Age Group and Virus Status")+
        theme_bw()+scale_fill_manual(values = phy_color)+theme(legend.position = "none")
p1+p2

### Season
df <- df4%>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
frequency_df <- df%>% group_by(season,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(season))%>%
        group_by(season, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(season, Virus, spe)  %>% 
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(season, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup()
ll <- frequency_df %>% filter(spe!="others") %>% group_by(spe,Virus) %>%  summarise(vv=mean(percentage)) %>% arrange(desc(vv)) %>% group_by(Virus) %>%
        mutate(rank = row_number(), 
               spe = ifelse(rank <= 15, spe, "others")) %>% filter(rank<16) %>% pull(spe)
frequency_df <- df%>% group_by(season,Virus,spe) %>% summarise(value=n()) %>% 
        filter(!is.na(season))%>%
        group_by(season, Virus, spe) %>%
        summarize(total_value = sum(value), .groups = 'drop') %>%
        arrange(desc(total_value)) %>%
        group_by(season, Virus, spe) %>% 
        summarize(total_value = sum(total_value), .groups = 'drop')%>%
        group_by(season, Virus) %>%
        mutate(percentage = total_value / sum(total_value) * 100) %>% ungroup() %>% filter(spe%in% ll) %>% mutate(spe=factor(spe,levels=rev(ll)))

## Southern northern detection diff
dd=df4   %>% filter(!duplicated(UID)) %>% select(south_north,type) %>% table %>% as.data.frame() 
df2 <- dd %>% group_by(south_north)%>%  mutate(csum = rev(cumsum(rev(Freq))), 
                                              pos = Freq/2 + lead(csum, 1),
                                              pos = if_else(is.na(pos), Freq/2, pos),
                                              value= (Freq/sum(Freq))*100)
df2$value <- round(df2$value,digits = 2)
 ggplot(df2, aes(x = "" , y = value, fill = fct_inorder(type))) +
  geom_col(width = 2, color = 1) +
  coord_polar(theta = "y") +facet_grid(.~south_north)+scale_fill_jama()+ theme_pubr()+xlab("")+ylab("")





## Co infection with HHV

ll=c("EBV","HCMV","HHV-7","HHV-1","HHV-3","HHV-6A")
df <- df4%>% filter(RBM=="Y",pathogen=="Y"|spe%in%ll) %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
frequency_df <- df%>% 
        group_by(UID) %>%filter(n()>1) %>%  do(data.frame(t(combn(.$spe %>% sort, 2)))) %>%
        ungroup() %>% dplyr::rename("Pathogen1"="X1", "Pathogen2"="X2") %>% dplyr::count(Pathogen1, Pathogen2) %>%
        right_join(table(df$spe) %>% as.data.frame() %>% dplyr::rename("Pathogen1"="Var1","total"="Freq"),.) %>% 
  dplyr::rename(Total1 = total) %>%
        right_join(table(df$spe) %>% as.data.frame() %>% dplyr::rename("Pathogen2"="Var1","total"="Freq"),.) %>% 
  dplyr::rename(Total2 = total) %>%
        mutate(Frequency = n / (Total1 + Total2)) %>%
        arrange(desc(Frequency)) %>% filter((Total2+Total1)>30) %>% filter(n>10) %>%
        filter(Pathogen1 %in% ll|Pathogen2 %in% ll) %>% filter(!(Pathogen1 %in% ll&Pathogen2 %in% ll))
for (i in 1:nrow(frequency_df)) {
        l1 <- frequency_df[i,"Pathogen2"]
        if (l1 %in% ll) {
                z=frequency_df[i,"Pathogen1"]
                frequency_df[i,"Pathogen2"]=z
                frequency_df[i,"Pathogen1"]=l1
        }
}

frequency_df$Frequency=frequency_df$n/frequency_df$Total2
wide_df <- frequency_df %>% dcast(Pathogen2~Pathogen1,value.var = "Frequency") %>% column_to_rownames("Pathogen2") %>% arrange(desc(HCMV))

# 使用dcast将长格式转换为宽格式，生成正方形矩阵
col_fun_prop = colorRamp2(c(0,0.05,0.1,0.15,0.2,0.3,0.4,0.5,0.6), 
                          colorRampPalette((RColorBrewer::brewer.pal(n =7, name = "YlGnBu")))(9))

ll1 <- rownames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus")) %>% arrange(group)
ll2 <- colnames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus")) %>% arrange(group)
Heatmap(wide_df[ll1$.,ll2$.],rect_gp = gpar(col="white",lwd=5),cluster_rows =F,cluster_columns = F,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = ll1$group,row_gap = unit(3, "mm"),column_split = ll2$group,
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,na_col = "grey87",
        right_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                          labels = c("Bacterial & Fungi","Virus"), 
                                                          labels_gp = gpar(col = "white", fontsize = 15))))

## hhv logistic
ll=c("EBV","HCMV","HHV-7","HHV-1","HHV-3","HHV-6A")
df <- df4%>% filter(RBM=="Y",pathogen=="Y"|spe%in%ll) %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
df <- df %>% dcast(season+Sexual+Age+Region~spe,fill = 0)
dim(df)
df1=df[,1:4]
df2=df[,5:60]
df2[df2!=0]=1
df=cbind(df1,df2)
ll=colnames(df) %>% .[-c(3)]
for (i in ll) {
        df[,i]=as.factor(df[,i])
}
library(marginaleffects)
library(gtsummary)
library(questionr)
ll=ll[-c(1:3)]
#test <- multinom(`Acinetobacter baumannii` ~ Sexual + Age+month + `Aspergillus flavus`, data = df)
dd=list()
z=1
### merge rsv flu 
ll2=c("EBV","HCMV","HHV-7","HHV-1","HHV-3","HHV-6A")
for (i in setdiff(ll,ll2)) {
        for (j in ll2) {
                test <- multinom(get(i) ~ Sexual + Age+Region+season + get(j), data = df)
                tt <- tidy(test, conf.int = TRUE) %>% as.data.frame()
                tt1 <- odds.ratio(test) %>% as.data.frame()
                tt$OR=tt1$OR
                tt$p=tt1$p
                tt$y=i
                tt$x=j
                dd[[z]]=tt
                z=z+1
                
        }
        
}
logistic_hsv <-  do.call(rbind,dd)

for (i in 1:nrow(frequency_df)) {
  j=frequency_df[i,"Pathogen2"]
    z=frequency_df[i,"Pathogen1"]
    dd[[i]]=logistic_hsv %>% filter(y==j,x==z) %>% filter(grepl("get",term)) %>% filter(!is.infinite(statistic))
}
df=do.call(rbind,dd) %>% mutate(p=p.adjust(p)) 

df=logistic_hsv %>% filter(grepl("get",term)) %>%mutate(p=p.adjust(p)) %>% filter(!is.infinite(statistic))

df=filter(df,p<0.05) 
wide_df <- df %>% filter(!duplicated(paste0(x,y))) %>% dcast(x~y,value.var = "OR") %>% column_to_rownames("x")
col_fun_prop = colorRamp2(c(0,0.05,0.1,0.3,1,1.7,1.9,3,4), 
                          colorRampPalette(rev(RColorBrewer::brewer.pal(n =8, name = "RdBu")))(9))

ll1 <- rownames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus")) %>% arrange(group)
ll2 <- colnames(wide_df) %>% as.data.frame() %>% mutate(group=ifelse(`.`%in% str_subset(ll_pathogen2,"virus",negate=T),"Bacterial","Virus")) %>% arrange(group)
Heatmap(wide_df[ll1$.,ll2$.],rect_gp = gpar(col="white",lwd=5),cluster_rows =F,cluster_columns = F,
        border = T,border_gp = gpar(col = "black", lwd=1.3),column_title_gp = gpar( col = "black", 
                                                                                    border = "black",lwd=2),
        row_split = ll1$group,row_gap = unit(3, "mm"),column_split = ll2$group,
        column_gap = unit(3, "mm"), 
        row_names_side = "left",
        col=col_fun_prop,na_col = "grey87",
        top_annotation = columnAnnotation(foo = anno_block(gp = gpar(fill = c("#374E55","#6B8636")),
                                                           labels = c("Bacterial & Fungi","Virus"), 
                                                           labels_gp = gpar(col = "white", fontsize = 15))))



ll=c("EBV","HCMV","HHV-7","HHV-1","HHV-3","HHV-6A")
df <- df4%>% filter(RBM=="Y",pathogen=="Y"|spe%in%ll) %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% filter((!is.na(PCT))|(!is.na(WBC)),(!is.na(CRP)))
#df=df %>% group_by(UID) %>% mutate(hsv=ifelse(length(unique(pathogen))==1,ifelse(all(spe %in% ll),"Only-HSV","Without-HSV"),"With-HSV")) %>% as.data.frame()

df=df %>% group_by(UID) %>% mutate(hsv=ifelse(length(unique(pathogen))==1,"Without-HHV","With-HHV")) %>% as.data.frame()
ll2=df %>% filter(hsv=="Only-HHV") %>% filter(!duplicated(UID)) %>% pull(UID)
## Top6
ll1=c("Pneumocystis jirovecii","Acinetobacter baumannii","Candida albicans","Rhinovirus A/B/C","RSV","HPIV-3","HCoV-229E")

ll=filter(df,spe=="Pneumocystis jirovecii") %>% pull(UID) %>% unique()
ggboxplot(df %>% filter(type!="No pathogen detected") %>% filter(UID %in% ll),"hsv","PCT",fill = "hsv")+stat_compare_means(label = "p.signif")
cc=list()
for (i in ll1) {
  ll=filter(df,spe==i) %>% pull(UID) %>% unique()
  cc[[i]] <- ggboxplot(df %>% filter(type!="No pathogen detected") %>% filter(UID %in% ll),
                       "hsv","PCT",color = "hsv",order = c("Without-HHV","With-HHV"),size=1.5)+
    stat_compare_means(method = "t.test",label = "p.signif")+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                                                                           labels = trans_format("log10", math_format(10^.x)))+
    theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))+ggtitle("",subtitle = i)+
    scale_color_aaas()+xlab("")#+geom_hline(yintercept = 0.25,color="grey45",size=1.5,linetype=2)
}
#filter(df,hsv=="Only-HSV") %>% filter(!is.na(PCT))%>% pull(PCT) %>% median()
p3 <- (purrr::reduce(cc, `+`) + plot_layout(nrow = 1))

p1/p2/p3



## DNA /RNA  ratio mNGS mTGS comparision
df=ABFV_s_de %>% rownames_to_column("spe") %>% right_join(abfv_list1,.) %>% filter(ABFV!="Virus") %>% select(-ABFV) %>% column_to_rownames("spe")
tt=df[,colSums(df)>300] %>% as.matrix() %>% prop.table(.,2) %>% as.data.frame()
df <- df4%>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
df$ratio=NA
df$dna=NA
df$rna=NA
for (i in 1:nrow(df)) {
        ll=df[i,"UID"]
        ll1 <- str_subset(ll_pair,paste0(ll,"-"))
        lr <- ll1[endsWith(ll1,"R")]
        ld <- ll1[endsWith(ll1,"D")]
        ss=df[i,"spe"]
        if (!is.null(tt[ss,lr])&(!is.null(tt[ss,ld]))) {
          df[i,"ratio"] <-  tt[ss,lr]/tt[ss,ld]
          df[i,"dna"] <-  tt[ss,ld]
          df[i,"rna"] <-  tt[ss,lr]
        }
}

dd1 <- df2_new %>% filter(spe %in% ll_new)
dd1 <- select(dd1,spe,UID)
dd1$dna=NA
dd1$rna=NA
for (i in 1:nrow(dd1)) {
  ll=dd1[i,"UID"]
  ll1 <- str_subset(ll_pair,paste0(ll,"-"))
  lr <- ll1[endsWith(ll1,"R")]
  ld <- ll1[endsWith(ll1,"D")]
  ss=dd1[i,"spe"]
  if (!is.null(tt[ss,lr])&(!is.null(tt[ss,ld]))) {
    dd1[i,"dna"] <-  tt[ss,ld]
    dd1[i,"rna"] <-  tt[ss,lr]
  }
}

dd2 <- select(df,spe   , UID     ,  dna    ,    rna) %>% rbind(.,dd1) %>% na.omit()
dd2$ratio=dd2$rna/dd2$dna
dd2$group=ifelse(dd2$spe %in% ll_new,"New","Old")
ll=dd2 %>% filter(dna!=0 & rna!=0)  %>% group_by(spe) %>% filter(n()>5) %>% mutate(ratio=log10(ratio)) %>%
  group_by(spe) %>% summarise(vv=median(ratio)) %>% arrange(vv) %>% pull(spe)
dd2$spe=factor(dd2$spe,levels = ll)
dd2$group=factor(dd2$group,levels = c("Old","New"))
dd2 %>% filter(dna!=0 & rna!=0) %>% filter(!is.na(spe)) %>% select(spe,dna,rna,group) %>% melt() %>%
  group_by(spe) %>% filter(n()>5) %>% mutate(value=log10(value+1e-4)) %>% 
  ggplot(.,aes(spe,value))+geom_boxplot(aes(fill=variable),alpha=0.8,outliers = F)+
  facet_grid(.~group,scales = "free_x",space = "free")+scale_fill_jama()+theme_pubr()+
  ylab("Pathogen relative abundance (Log10)")+
  stat_compare_means(aes(group=variable),label = "p.signif")+
  theme(axis.text.x = element_text(angle =45,hjust = 1))
  
  ggboxplot(.,order = ll, x = "spe",y="value",fill = "variable",alpha=0.8,outlier.shape = NA,facet.by = "group")+
  coord_flip()+scale_fill_jama()+ylab("Pathogen relative abundance (Log10)")+
  stat_compare_means(aes(group=variable),label = "p.signif")

df %>% filter(dna!=0 & rna!=0) %>% select(spe,dna,rna) %>% melt() %>% group_by(spe) %>% filter(n()>15) %>% mutate(value=log10(value+1e-5)) %>% 
  ggboxplot(.,order = ll, x = "spe",y="value",fill = "variable",alpha=0.8,outlier.shape = NA)+
  coord_flip()+scale_fill_jama()+ylab("Pathogen relative abundance (Log10)")+
  stat_compare_means(aes(group=variable),label = "p.signif")

# all  at least greater than 1% in one sample
dd=list()
for (i in gsub("-.*","",ll_pair) %>% unique()) {
        ll1 <- str_subset(ll_pair,paste0(i,"-"))
        lr <- ll1[endsWith(ll1,"R")]
        ld <- ll1[endsWith(ll1,"D")]
        if (all(c(ld,lr) %in% colnames(tt))) {
                a <- tt[,c(lr,ld)] %>% .[rowSums(.>0.01)>0,] %>% .[rowSums(.> 1e-5)>1,] %>% rownames_to_column("spe") %>% mutate(ratio=get(lr)/get(ld),diff=(get(lr)-get(ld)),UID=i) %>% na.omit() %>% filter(!is.infinite(ratio))  %>% mutate(rank_r=rank(ratio),rank_d=rank(diff))
                colnames(a)[2:3]=c("RNA","DNA")
                dd[[i]] <- a
        }
}
ratio_all <-  do.call(rbind,dd)

ratio1 <- select(df4,UID,spe,RBM,pathogen) %>% right_join(.,ratio_all) 
ratio1$pathogen[is.na(ratio1$pathogen)]="N"

tt1=ratio1 %>% filter(pathogen=="Y") %>% group_by(spe) %>% filter(n()>15) %>% mutate(DNA=log10(DNA+1e-5),RNA=log10(RNA+1e-5)) %>% arrange(spe) %>% select(spe,DNA,RNA) 
ll=ratio1 %>% filter(pathogen=="Y") %>% group_by(spe) %>% filter(n()>15) %>% mutate(ratio=log10(ratio)) %>% group_by(spe) %>% summarise(vv=median(ratio)) %>% arrange(vv) %>% pull(spe)
  ggboxplot(tt1,order = ll,
            x = "spe",y="value",fill = "variable",alpha=0.8)+
  coord_flip()+scale_fill_jama()+ylab("Pathogen relative abundance (Log10)")+
    stat_compare_means(aes(group=variable),label = "p.signif")+ylim(c(-3,1))




library(FFTrees)
ratio2 <- filter(ratio1,spe%in% ll_pathogen2) %>% right_join(metadata %>% filter(!duplicated(UID)),.)
table(ratio2$spe) %>% sort %>% tail(20)
df=filter(ratio2,spe=="Mycoplasmoides pneumoniae")

df1 <- df %>% select(pathogen,Department,Sexual,Age,month,ratio,diff,rank_r,rank_d) %>% na.omit()
df1$pathogen=as.factor(df1$pathogen)
df1$Department=as.factor(df1$Department)
df1$Sexual=as.factor(df1$Sexual)
#df1$Province=as.factor(df1$Province)
trn <- sample(1:nrow(df1), size = floor(nrow(df1)*0.75), replace = FALSE)
heart_fft <- FFTrees(formula = pathogen ~., algorithm = "dfan",
                     data = df1[trn,],
                     data.test = df1[-trn,], 
                     decision.labels = c( "No","Pathogen"))
plot(heart_fft,data = "test")

###  Ratio 
ll1 <- table(ratio1$spe) %>% sort() %>% tail(112) %>% as.data.frame() %>%
        mutate(group=ifelse(Var1 %in% ll_pathogen2,"Pathogen","No")) %>% arrange(group,desc(Freq)) %>%
        group_by(group) %>% mutate(rank=1:n()) %>%
        filter(rank<22) %>% filter(group=="Pathogen") %>% as.data.frame() %>% pull(Var1) %>% as.character()
ll2 <- table(ratio1$spe) %>% sort() %>% tail(112) %>% as.data.frame() %>%
        mutate(group=ifelse(Var1 %in% ll_pathogen2,"Pathogen","No")) %>% arrange(group,desc(Freq)) %>%
        group_by(group) %>% mutate(rank=1:n()) %>%
        filter(rank<22) %>% filter(group=="No") %>% as.data.frame() %>% pull(Var1) %>% as.character()
ll1 <- ratio1 %>% filter(spe %in% ll1)%>% mutate(ratio=log10(ratio) ) %>% group_by(spe) %>% summarise(aa=median(ratio)) %>% arrange(aa) %>% pull(spe)
ll2 <- ratio1 %>% filter(spe %in% ll2)%>% mutate(ratio=log10(ratio) ) %>% group_by(spe) %>% summarise(aa=median(ratio)) %>% arrange(aa) %>% pull(spe)
p1 <- ggboxplot(ratio1 %>% filter(spe %in% ll1)%>% mutate(ratio=log10(ratio)) %>% arrange(spe),
          x = "spe",y="ratio",fill = "pathogen",order = ll1,alpha=0.8)+
        stat_compare_means(aes(group=pathogen),label = "p.signif")+
        coord_flip()+scale_fill_jama()+ylab("Log10( RNA/DNA )")+
        geom_hline(yintercept = 0,linetype=2)

p2 <- ggboxplot(ratio1 %>% filter(spe %in% ll2)%>% mutate(ratio=log10(ratio)) %>% arrange(spe),
          x = "spe",y="ratio",fill = "#374E55",order = ll2,alpha=0.8)+
        coord_flip()+scale_fill_jama()+ylab("Log10( RNA/DNA )")+
        geom_hline(yintercept = 0,linetype=2)

p1+p2

##new

df=ratio_all %>% filter(spe %in% ll_pathogen2)
ll1 <- df %>% filter(spe %in% ll_pathogen2)%>% mutate(ratio=log10(ratio) ) %>% group_by(spe) %>% summarise(aa=median(ratio)) %>% arrange(aa) %>% pull(spe)
df %>% mutate(ratio=log10(ratio+1e-5)) %>% ggboxplot(.,x = "spe",y="ratio",fill = "#374E55",order = ll1,alpha=0.8)+coord_flip()+scale_fill_jama()+ylab("Log10( RNA/DNA )")+
  geom_hline(yintercept = 0,linetype=2)




#### RF
ll <- table(ratio2$spe,ratio2$pathogen) %>% as.data.frame() %>% dcast(Var1~Var2) %>% arrange(desc(Y)) %>% head(21) %>% 
        filter(Var1!="Mycobacterium tuberculosis") %>% filter(N >10) %>% pull(Var1) %>% as.character()
ll1 <- filter(ratio2,spe%in%ll)  %>% compare_means(data = .,ratio~pathogen,group.by = "spe") %>% filter(p.adj<0.05) %>% pull(spe)
ggboxplot(filter(ratio2,spe%in%ll1) %>% mutate(ratio=log10(ratio)),x = "spe",y="ratio",fill = "pathogen")+
        stat_compare_means(aes(group=pathogen),label = "p.signif")+coord_flip()+scale_fill_jama()+ylab("Log10( RNA/DNA )")+geom_hline(yintercept = 0,linetype=2)

ll1 <- filter(ratio2,spe%in%ll)  %>% compare_means(data = .,diff~pathogen,group.by = "spe") %>% filter(p.adj<0.05) %>% pull(spe)
ggboxplot(filter(ratio2,spe%in%ll1) %>% mutate(ratio=log10(ratio)),x = "spe",y="diff",fill = "pathogen")+
        stat_compare_means(aes(group=pathogen),label = "p.signif")+coord_flip()+scale_fill_jama()
df=filter(df4,type=="No pathogen detected")  %>% filter((krakenD>2&krakenR>2)) %>% filter(spe!="Pseudomonas aeruginosa",spe!="Stenotrophomonas maltophilia")
ll=table(df$spe) %>% sort %>% as.data.frame() %>% tail(20) %>% filter(!(Var1 %in% c("HHV-1","EBV","HCMV"))) %>% pull(Var1) %>% as.character()
filter(ratio1,spe%in% ll) %>% right_join(metadata %>% filter(!duplicated(UID)),.) %>% mutate(ratio=log10(ratio))%>% 
         filter(grepl("Prevo",spe)|grepl("Veill",spe)|grepl("Rothia",spe)) %>% arrange(spe) %>% 
        ggboxplot(.,x = "spe",y="ratio",fill = "#374E55")+coord_flip()+ylab("Log10( RNA/DNA )")+geom_hline(yintercept = 0,linetype=2)


### pathogen rank
df=ABFV_s_p
df[df<0.005]=0
df1 <- apply(df, 2, function(x) rank(-x, ties.method = "min"))  
ll=intersect(ll_pathogen2,rownames(df1))
df1 <- df1[ll,]
replace_most_frequent_with_zero <- function(column) {
  # 计算频率
  freq <- table(column)
  # 找出最频繁的数字
  most_frequent <- as.numeric(names(freq)[which.max(freq)])
  # 将最频繁的数字替换为0
  column[column == most_frequent] <- 0
  return(column)
}

# 应用这个函数到每一列
df1=matrix_data <- apply(df1, 2, replace_most_frequent_with_zero)
df2=df1 %>% as.data.frame() %>% rownames_to_column("spe") %>% melt %>% filter(value!=0) %>% group_by(variable) %>% 
  filter(n()>=2) %>% arrange(variable,value) %>% mutate(rank=1:n())
df2=right_join(abfv_list1,df2) %>% mutate(ABFV=gsub("Fungi","Bacterial",ABFV)) %>% group_by(variable,ABFV) %>%filter(n()>=2) %>%  mutate(rank=1:n())
species_average_rank <- df2 %>%
  group_by(spe,ABFV) %>%
  summarise(AverageRank = mean(rank, na.rm = TRUE)) %>% filter(AverageRank<2)%>% arrange(ABFV,AverageRank) 
species_average_rank$spe <- gsub("Human gammaherpesvirus 4","EBV",species_average_rank$spe)
species_average_rank$spe <- gsub("Human betaherpesvirus 5","HCMV",species_average_rank$spe)
species_average_rank$spe <- gsub("Human betaherpesvirus 7","HHV-7",species_average_rank$spe)
species_average_rank$spe <- gsub("Human alphaherpesvirus 1","HHV-1",species_average_rank$spe)
species_average_rank$spe <- gsub("Human alphaherpesvirus 3","HHV-3",species_average_rank$spe)
species_average_rank$spe <- gsub("Human betaherpesvirus 6A","HHV-6A",species_average_rank$spe)
species_average_rank$spe <- gsub("Rhinovirus .*","Rhinovirus A/B/C",species_average_rank$spe)
species_average_rank$spe <- gsub("Orthopneumovirus hominis","RSV",species_average_rank$spe)
species_average_rank$spe <- gsub("Respirovirus laryngotracheitidis","HPIV-1",species_average_rank$spe)
species_average_rank$spe <- gsub("Influenza .* virus","Influenza A/B/C",species_average_rank$spe)
species_average_rank$spe <- gsub("Respirovirus pneumoniae","HPIV-3",species_average_rank$spe)
species_average_rank$spe <- gsub("Human orthorubulavirus 4","HPIV-4",species_average_rank$spe)
species_average_rank$spe <- gsub("Severe acute respiratory syndrome-related coronavirus","SARS-CoV-2",species_average_rank$spe)
species_average_rank$spe <- gsub("Enterovirus .*","Enterovirus",species_average_rank$spe)
species_average_rank$spe <- gsub("Human coronavirus 229E","HCoV-229E",species_average_rank$spe)
species_average_rank$spe <- gsub("Betacoronavirus 1","HCoV-OC43",species_average_rank$spe)
species_average_rank$spe <- gsub("Human coronavirus HKU1","HCoV-HKU1",species_average_rank$spe)
species_average_rank$spe <- gsub("Human coronavirus NL63","HCoV-NL63",species_average_rank$spe)
species_average_rank$spe <- gsub("Human mastadenovirus .*","HAdV",species_average_rank$spe)
species_average_rank$spe <- gsub("Metapneumovirus hominis","HMPV",species_average_rank$spe)
species_average_rank <- species_average_rank %>% group_by(spe,ABFV) %>% summarise(AverageRank=mean(AverageRank)) %>% as.data.frame() %>% 
  arrange(ABFV,AverageRank)
###
df=filter(df4,spe=="Influenza A/B/C"|spe=="Staphylococcus aureus") %>%  filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2)))
ll=df$UID[duplicated(df$UID)]
df=filter(df4,spe=="Staphylococcus aureus") %>%  filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% mutate(group=ifelse(UID %in% ll,"Co-HCMV","Fungi"))
ggboxplot(df,"group","CRP")+stat_compare_means()


###new pathogen
ll=df4 %>% filter(type=="No pathogen detected") %>% pull(UID) %>% 
  c(.,ll_pair %>% as.data.frame() %>% dplyr::rename("ID"=".") %>%
  mutate(UID=gsub("-.*","",ID)) %>% filter(!(UID %in% df4$UID)) %>% pull(UID)) %>% 
  unique()
df=ABFV_s_de %>% rownames_to_column("spe") %>% right_join(abfv_list1,.) %>% filter(ABFV!="Virus") %>% dplyr::select(-ABFV) %>% column_to_rownames("spe")
tt=df[,colSums(df)>200] %>% as.matrix() %>% prop.table(.,2) %>% as.data.frame()

df=tt %>% rownames_to_column("spe") %>% melt %>% filter(value>0.4) %>% 
  filter(!(spe %in% ll_pathogen2)) %>% mutate(UID=gsub("-.*","",variable)) %>% filter(UID %in% ll) %>% filter(variable %in% ll_pair)
p1 <- table(df %>% filter(!duplicated(UID)) %>% pull(spe)) %>% as.data.frame() %>% filter(Freq>8) %>% arrange(desc(Freq))%>% ggbarplot(.,"Var1","Freq",sort.val = "asc",fill = "#0E2A59")+coord_flip()
df1=filter(df,spe %in% as.character(p1$data$Var1)) %>%filter(variable %in% ll_pair) %>%  distinct(spe,UID,.keep_all = T) %>% right_join(metadata1,.)
df1$paired_value=NA
for (i in 1:nrow(df1)) {
  ll2 <- str_subset(ll_pair,paste0(df1[i,"UID"],"-"))
  lr <- setdiff(ll2,df1[i,"variable"])
  if (lr %in% colnames(tt)) {
    df1[i,"paired_value"]=tt[df1[i,"spe"],lr]
  }
}
df1=filter(df1,paired_value>0)

ll=df1 %>% filter(!duplicated(paste0(UID,spe))) %>% pull(UID)
tt1=filter(metadata %>% filter(ID %in% ll_pair),UID %in% ll) 
tt2=ABFV_s_p_de[ll_pathogen2,tt1$ID] %>% rownames_to_column("spe") %>% melt %>%
  filter(value>0.01) %>% filter(!grepl("herpesvirus",spe)) %>% arrange(desc(value))
tt2$paired_value=NA
tt2$UID=gsub("-.*","",tt2$variable)
for (i in 1:nrow(tt2)) {
  ll2 <- str_subset(ll_pair,paste0(tt2[i,"UID"],"-"))
  lr <- setdiff(ll2,tt2[i,"variable"])
    tt2[i,"paired_value"]=ABFV_s_p_de[tt2[i,"spe"],lr]
}
ll_newpathogen_co_id <-  tt2 %>% filter(value>0.05,paired_value>0.05) %>% mutate(UID=gsub("-.*","",variable)) %>% filter(!duplicated(UID)) %>% pull(UID)



df_new_pathogen=filter(df1 %>% filter(!duplicated(UID)),paired_value>0) %>% filter(!(UID %in% ll_newpathogen_co_id)) %>% filter(spe%in%ll_new)

df2=df1 %>% filter(!duplicated(paste0(UID,spe))) %>%  filter(!(UID %in% ll_newpathogen_co_id))

ll3=table(df2$spe) %>% as.data.frame() %>% arrange(Freq) %>% filter(Freq>=7)
l1=df2 %>% group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>3) %>% filter(!is.na(WBC)) %>% filter(!(UID %in% ll_newpathogen_co_id))%>% summarise(tt=median(WBC)) %>% arrange(desc(tt)) %>% pull(spe) %>% head(15)
l2 <- df2 %>% group_by(spe) %>% mutate(nn=sum(!is.na(CRP))) %>% filter(nn>3) %>% filter(!is.na(CRP)) %>% filter(!(UID %in% ll_newpathogen_co_id))%>% summarise(tt=median(CRP)) %>% arrange(desc(tt)) %>% pull(spe) %>% head(15)
l3 <- df2 %>% group_by(spe) %>% mutate(nn=sum(!is.na(PCT))) %>% filter(nn>3) %>% filter(!is.na(PCT)) %>% filter(!(UID %in% ll_newpathogen_co_id))%>% summarise(tt=median(PCT)) %>% arrange(desc(tt)) %>% pull(spe) %>% head(15)
c(l1,l2,l3) %>% table %>% sort %>% as.data.frame() %>% filter(`.` %in% ll3$Var1)
ll_new=c(l1,l2,l3) %>% table %>% sort %>% as.data.frame() %>% filter(`.` %in% ll3$Var1)%>% filter(Freq>1) %>% pull(1) %>% as.character()
ll_new=ll_new[which(!(ll_new %in% c("Sphingobium sp. YG1","Streptococcus oralis","Rothia mucilaginosa","Ralstonia pickettii")))]
ll_new
setdiff(ll_new_old,ll_new)
setdiff(ll_new,ll_new_old)
dd=filter(df2,spe %in% ll_new)
dd$other=0
dd$other1=NA

ll1=ll_pathogen2[!grepl("virus",ll_pathogen2)]
df=ABFV_s_p_de[ll1,]
for (i in 1:nrow(dd)) {
  ll2 <- str_subset(ll_pair,paste0(dd[i,"UID"],"-"))
  dd[i,"other"]=df %>% na.omit() %>% select(ll2) %>% .[rowSums(.>0.01)>1,] %>% nrow()
  dd[i,"other1"]=df %>% na.omit() %>% select(ll2) %>% .[rowSums(.>0.01)>1,] %>% rownames() %>% paste(.,collapse = "_")
  
}
psych::corr.test(ABFV_s_p_de["Candida albicans",] %>% select(ends_with("R")) %>% as.numeric(),ABFV_s_p_de["Nakaseomyces glabratus",] %>% select(ends_with("R")) %>% as.numeric(),method = "spearman")
ll_newpathogen_co_id <- filter(dd,other!=0,spe!="Nakaseomyces glabratus") %>% pull(UID) %>% unique()


tt1 <- table(df1 %>% filter(!duplicated(paste0(UID,spe))) %>% pull(spe)) %>% as.data.frame() %>% filter(Freq>10) %>% arrange((Freq)) 
p1 <- ggbarplot(tt1,"Var1","Freq",sort.val = "asc",fill = "#385581")+coord_flip()+ theme(axis.text.y = element_text(colour = ifelse(tt1$Var1 %in% ll_new,"red","black")))
ll1 <- df2 %>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa")))%>% group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>4) %>% filter(!is.na(WBC)) %>% summarise(tt=median(WBC)) %>% arrange((tt)) %>% pull(spe)
#ggboxplot(df1 %>% group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>5) %>% mutate(tt=median(value)) %>% arrange(desc(tt)),"spe","WBC")+coord_flip()
p2 <- ggboxplot(df2 %>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa")))%>%  group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>4) %>% filter(!is.na(WBC)) %>% filter(!(UID %in% ll_newpathogen_co_id)),"spe","WBC",order = ll1,fill = "#dac190")+geom_hline(aes(yintercept=10),colour="grey55", linetype="dashed",size=1)+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+coord_flip()+
  theme(axis.text.y = element_text(colour = ifelse(ll1 %in% ll_new,"red","black")))
ll1 <- df2 %>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa")))%>%
  group_by(spe) %>% mutate(nn=sum(!is.na(CRP))) %>% filter(nn>3) %>% filter(!is.na(CRP)) %>% summarise(tt=median(CRP)) %>% arrange((tt)) %>% pull(spe)
#ggboxplot(df1 %>% group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>5) %>% mutate(tt=median(value)) %>% arrange(desc(tt)),"spe","WBC")+coord_flip()
p3 <- ggboxplot(df2  %>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa")))%>% 
                  group_by(spe) %>% mutate(nn=sum(!is.na(CRP))) %>% filter(nn>3) %>% filter(!is.na(CRP))%>% filter(!(UID %in% ll_newpathogen_co_id)),"spe","CRP",order = ll1,fill = "#dac190")+geom_hline(aes(yintercept=10),colour="grey55", linetype="dashed",size=1)+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),   labels = trans_format("log10", math_format(10^.x)))+coord_flip()+
  theme(axis.text.y = element_text(colour = ifelse(ll1 %in% ll_new,"red","black")))

ll1 <- df2 %>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa"))) %>%
  group_by(spe) %>% mutate(nn=sum(!is.na(PCT))) %>% filter(nn>3) %>% filter(!is.na(PCT)) %>% summarise(tt=median(PCT)) %>% arrange((tt)) %>% pull(spe)
#ggboxplot(df1 %>% group_by(spe) %>% mutate(nn=sum(!is.na(WBC))) %>% filter(nn>5) %>% mutate(tt=median(value)) %>% arrange(desc(tt)),"spe","WBC")+coord_flip()
p4 <- ggboxplot(df2%>% filter(!(spe %in% c("Streptococcus oralis","Rothia mucilaginosa"))) %>%
                  group_by(spe) %>% mutate(nn=sum(!is.na(PCT))) %>% filter(nn>3)%>% filter(!(UID %in% ll_newpathogen_co_id)),"spe","PCT",order = ll1,fill = "#dac190")+geom_hline(aes(yintercept=0.5),colour="grey55", linetype="dashed",size=1)+scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),  labels = trans_format("log10", math_format(10^.x)))+coord_flip()+
  theme(axis.text.y = element_text(colour = ifelse(ll1 %in% ll_new,"red","black")))




p1|p2|p3|p4
ll2=c("Candida parapsilosis","Nakaseomyces glabratus","Metamycoplasma hominis","Carnobacterium maltaromaticum","Corynebacterium striatum","Parvimonas micra","Acinetobacter junii","Tropheryma whipplei")


df1=ABFV_s_p_de[c(ll_pathogen2,ll_new),] %>% select(ends_with("D")) %>% na.omit()
df2=ABFV_s_p_de[c(ll_pathogen2,ll_new),] %>% select(ends_with("R")) %>% na.omit()
library(psych)
tt=metadata %>% filter(!is.na(PCT),ID %in% colnames(df1))
tt1=df1[,tt$ID]
tt2=corr.test(tt[,"WBC"],t(tt1),method = "spearman")
tt2$r %>% t %>% as.data.frame() %>% arrange(V1)



## VF stastical
vf= read_table("Downloads/all_vf1.tsv",col_names = F)
colnames(vf)=c("ID","VFGID","length","counts")
vf1 <- vf %>% filter(counts*25>length)%>% dcast(VFGID+length~ID)
vf_ann <- read_csv("Downloads/vf_ann2.csv",col_names = T)
vf1[is.na(vf1)] <- 0
vf2 <- dplyr::select(vf1,-2) %>% right_join(dplyr::select(vf_ann,1,3),.)
vf1[1:5,1:5]


ll=ll_pathogen2[!grepl("virus",ll_pathogen2)]
df=tt[ll,intersect(ll_pair,colnames(tt))] %>% na.omit %>% rownames_to_column("spe") %>% melt
colnames(df)[2]="ID"
df$UID=gsub("-.*","",df$ID)
df=filter(df,value!=0)
dd=df %>% group_by(UID,spe) %>% summarise(value=max(value)) 
dd$pathogen="N"
dd=as.data.frame(dd)
df <- df4%>% filter(RBM=="Y",pathogen=="Y") %>% filter((Virus=="N"&(gapD>1|gapR>1)&(krakenD>2&krakenR>2))|(Virus=="Y"&(gapD>1|gapR>1)&(krakenD>2|krakenR>2))) %>% filter(Virus=="N")
for (i in 1:nrow(df)) {
  dd[which((dd$UID==df[i,"UID"])&(dd$spe==df[i,"spe"])),"pathogen"]="Y"
}
dd1 <- dd %>% group_by(spe,pathogen) %>% summarise(Q1 = quantile(value, 0.25),
                                            Q3 = quantile(value, 0.75),
                                            IQR = median(value),n=n()) %>%
  pivot_wider(
    names_from = pathogen,
    values_from = c(n,Q1, Q3, IQR),
    names_sep = "_"
  ) 
write.csv(dd1,"pathogen_abundance.csv")


df=ABFV_s_de %>% rownames_to_column("spe") %>% right_join(abfv_list1,.) %>% filter(ABFV!="Virus") %>% select(-ABFV) %>% column_to_rownames("spe")
tt=df[,colSums(df)>500] %>% as.matrix() %>% prop.table(.,2) %>% as.data.frame()