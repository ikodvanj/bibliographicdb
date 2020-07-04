library(dplyr)
library(ggplot2)
library(stringr)
library(devEMF)
library(ggsci)
library(pubmedR)
library(bibliometrix)
library(naniar)
library(VennDiagram)

theme_set(theme_bw() + 
            theme( plot.subtitle = element_text(vjust = 1), 
                   plot.caption = element_text(vjust = 1), 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(), 
                   axis.title = element_text(size = 9), 
                   axis.text = element_text(size = 8, colour = "black"), 
                   panel.background = element_rect(fill = NA),
                   panel.border = element_rect(linetype = "solid", fill = NA, size = 1))
)

deduplicate <- function(s){
  return(s[!duplicated(s)])
}

Venn <- function(who, pubmed, scopus, searchphrase){
  df <- who
  df$Fulltext.URL[is.na(df$Fulltext.URL) & !is.na(df$DOI)] <- df$DOI[is.na(df$Fulltext.URL) & !is.na(df$DOI)]
  df$Fulltext.URL <- gsub("https://doi.org/", "", df$Fulltext.URL)
  df$Fulltext.URL <- gsub("https://dx.doi.org/", "", df$Fulltext.URL)
  df$Fulltext.URL[df$Fulltext.URL == ""] <- NA
  venn.diagram(
    x = list(deduplicate(na.omit(df$Fulltext.URL)), 
             deduplicate(na.omit(pubmed$DI)), 
             deduplicate(na.omit(scopus$DOI))),
    category.names = c("WHO" , "PubMed " , "Scopus"),
    col=c("#00B4E2FF", '#5C88DAFF', '#FFCD00FF'),
    fill = c(alpha("#00B4E2FF",0.3), alpha('#5C88DAFF',0.3), alpha('#FFCD00FF',0.3)),
    imagetype="png",
    # height = 1600, 
    # width = 1600,
    resolution = 900,
    filename = paste("venn_",searchphrase, ".png", sep = ""),
    output=TRUE
    
    )
}

updateData <- function(who, pubmed, scopus, topic, data){
  return(
    rbind(data,
          data.frame(Topic = topic, who = nrow(who), pubmed = nrow(pubmed), scopus = nrow(scopus)))
  )
}

clean <- function(df){
  df[df=="NA"]<-NA
  df[df==""]<-NA
  df[df==" "]<-NA
  df[df=="  "]<-NA
  df[df=="   "]<-NA
  df[df=="    "]<-NA
  df[df=="N/A"]<-NA
  df[df=="N/A."]<-NA
  df[df=="None"]<-NA
  df[df=="None."]<-NA
  df[df=="none"]<-NA
  df[df=="Not available"]<-NA
  df[df=="Not available."]<-NA
  df[df=="[No abstract available]"]<-NA
  df[df=="No abstract available."]<-NA
  return(df)
}

pubmedformating <- function(topic){
  search <- paste(topic, " AND COVID-19", sep = "")
  pubmed <- pmApiRequest(search, limit=21000, api_key = "") # INSERT API KEY
  pubmed <- pmApi2df(pubmed, format = "bibliometrix")
  
  pubmed <- pubmed %>% mutate_all(as.character())
  pubmed <- clean(pubmed)
  
  return(pubmed)
}

miss <- function(df, df2, fill, filter, names){
  dimensions <- dim(df)
  p1 <- miss_case_summary(df) %>%
    group_by(n_miss) %>%
    summarise(count = n()) %>%
    mutate(percentage = count / dimensions[1] * 100,
           n_miss_percentage = round(n_miss / dimensions[2] * 100, 0),
           x = paste(n_miss, "\n" , n_miss_percentage, sep = "")) %>%
    ggplot(aes(x = reorder(factor(x), n_miss), y = percentage, fill = "Color")) + geom_col() + 
      xlab("No. & % of NA per entry [n, %]") + ylab("Percentage of entries [%]")+
      scale_fill_manual(values = c(fill)) + theme(legend.position = "None") + 
      scale_y_continuous(limits = c(0,35)) 
  p2 <- miss_var_summary(df) %>%
    ggplot(aes(x = reorder(variable, pct_miss), y = pct_miss, fill = "Color")) + geom_col() + scale_y_continuous(limits = c(0,100)) + 
      xlab("") + ylab("Percentage of missing values [%]") + coord_flip()+
      scale_fill_manual(values = c(fill)) + theme(legend.position = "None")
  
  p4 <- vis_miss(df) + theme(axis.text.x = element_text(angle = 0)) + coord_flip()
  
  y <- list()
  
  for (i in 1:dim(df2)[2]){
    y[i] <- sum(duplicated(na.omit(df2[,i])) | duplicated(na.omit(df2[,i]), fromLast=TRUE)) / length(na.omit(df2[,i])) * 100
  }
  dfDup<-data.frame(dup = unlist(y))
  
  dfDup$category <- colnames(df2)[1:dim(df2)[2]]
  
  p3 <- dfDup[dfDup$category %in% filter,] %>%
    ggplot( aes(x = category, y = dup, fill = "color")) + geom_col() + 
    scale_y_continuous(limits = c(0,25)) + 
    xlab("") + ylab("Percentage duplicated data [%]") +
    scale_x_discrete(labels= names) +
    scale_fill_manual(values = c(fill)) + theme(legend.position = "None") 
  
  return(list(p1 = p1, p2 = p2, p3 = p3, p4 = p4))
}


miss_spec <- function(df){

  p4 <- vis_miss(df, cluster = TRUE) + theme(axis.text.x = element_text(angle = 0)) + coord_flip()

  return(p4)
}

extract = function(){}

matching = function(who, scopus, pubmed, topic){
  
  p = select(pubmed, AB, TI, ID,DE, DI)
  w = select(who, Abstract, Title, Descriptor.s., DOI, Fulltext.URL)
  s = select(scopus, Abstract, Title, Index.Keywords, Author.Keywords, DOI)
  
  for(i in 1:4){
    p[,5+i] = str_detect(toupper(p[,i]), toupper(topic))
  }
  pf = select(p, -AB, -TI, -ID, -DE)
  colnames(pf) = c("DI", "AB", "TI", "ID", "DE")
  
  for(i in 1:3){
    w[,5+i] = str_detect(toupper(w[,i]), toupper(topic))
  }
  wf = select(w, -Abstract, -Title, -Descriptor.s.,)
  colnames(wf) = c("DI", "URL", "Abstract", "Title", "Descriptors")
  
  for(i in 1:4){
    s[,5+i] = str_detect(toupper(s[,i]), toupper(topic))
  }
  sf = select(s, -Abstract, - Title, - Index.Keywords, - Author.Keywords)
  colnames(sf) = c("DI", "AB", "TI", "ID", "DE")
  
  pf = mutate(pf, topic = topic)
  sf = mutate(sf, topic = topic)
  wf = mutate(wf, topic = topic)
  
  return(list(pf,sf,wf))
  
}

plotting2 = function(l, topic){
  pmx =  l[[1]] 
  sfx =  l[[2]]
  wfx =  l[[3]]

  pmx[is.na(pmx)]<-FALSE
  sfx[is.na(sfx)]<-FALSE
  wfx[is.na(wfx)]<-FALSE
  
  pmx[pmx == TRUE]<-NA
  sfx[sfx == TRUE]<-NA
  wfx[wfx == TRUE]<-NA
  
  pmx = pmx %>% select(AB, TI, ID, DE)
  sfx = sfx %>% select(AB, TI, ID, DE)
  wfx = wfx %>% select(Abstract, Title, Descriptors)
  
  return(list(
    p1 = miss_spec(pmx) + scale_fill_manual(values = c("#F7F7F7", "#5C88DAFF") ,labels = c("Missmatch", "Match")) + theme(legend.title = element_blank()),
    p2 = miss_spec(sfx) + scale_fill_manual(values = c("#F7F7F7", "#FFCD00FF") ,labels = c("Missmatch", "Match")) + theme(legend.title = element_blank()),
    p3 = miss_spec(wfx) + scale_fill_manual(values = c("#F7F7F7", "#00B4E2FF") ,labels = c("Missmatch", "Match")) + theme(legend.title = element_blank())
  ))
}

saveplot = function(xplot, topicx){
  emf(paste("st_", topicx, "_pubmed.emf", sep = ""), width = 3.85, height = 2.4)
  print(xplot$p1)
  dev.off()
  emf(paste("st_", topicx, "_scopus.emf", sep = ""), width = 3.85, height = 2.4)
  print(xplot$p2)
  dev.off()
  emf(paste("st_", topicx, "_who.emf", sep = ""), width = 3.85, height = 2.4)
  print(xplot$p3)
  dev.off()
  
  tiff(paste("st_", topicx, "_pubmedL.tiff", sep = ""), width = 3.85*400, height = 2.4*400)
  print(xplot$p1)
  dev.off()
  tiff(paste("st_", topicx, "_scopusL.tiff", sep = ""), width = 3.85*400, height = 2.4*400)
  print(xplot$p2)
  dev.off()
  tiff(paste("st_", topicx, "_whoL.tiff", sep = ""), width = 3.85*400, height = 2.4*400)
  print(xplot$p3)
  dev.off()
}





df <- readRDS("WHOdb.RDS")

df <- clean(df)

cplot1 <- miss(
  df,
  df %>% mutate(Title.Authors.Journal = paste(Title, Authors, Journal)),
  fill = "#00B5E2FF", 
  filter = c("Abstract", "DOI", "Fulltext.URL", "ID", "Title.Authors.Journal"),
  names = c("Abstract", "DOI", "FullText URL", "ID", "Title\nAuthors\nJournal")
)

emf("who1.emf", width = 3, height = 2.5)
print(cplot1$p1)
dev.off()

emf("who2.emf", width = 3, height = 3.7)
print(cplot1$p2)
dev.off()

emf("who3.emf", width = 3, height = 2.7)
print(cplot1$p3)
dev.off()

emf("who4.emf", width = 5.8, height = 2.9)
print(cplot1$p4)
dev.off()

tiff("who4t.tiff", width = 5.8 * 413, height = 2.9 * 413)
print(cplot1$p4)
dev.off()


df$ID <- seq(1, dim(df)[1], 1)
df_url <- df[is.na(df$Fulltext.URL) != TRUE,] # excludamo duplikate po DOI-u tako da kasnije mozemo traziti duplikate po DOI-u
table(duplicated(df_url$Fulltext.URL) | duplicated(df_url$Fulltext.URL, fromLast=TRUE))
duplicatedurl <- df_url[duplicated(df_url$Fulltext.URL) | duplicated(df_url$Fulltext.URL, fromLast=TRUE),]
final1 <- rbind(duplicatedurl)
final1 <- final1[duplicated(final1$ID) == FALSE,]
final1 %>%
  select(ID, Authors, Title, "Journal", DOI, Fulltext.URL, Abstract) %>%
  arrange(Fulltext.URL) %>%
  xlsx::write.xlsx("WHO - final db - duplicated based on URL.xlsx", row.names = FALSE)









## Scopus ###########################################################################
## Used search phrase:                                                              #
## "COVID-19"                                                                       #  
##                                                                                  #
## Accessed on: 2020-06-26 ##########################################################


# scopus <- convert2df("scopus_top2000.bib", dbsource = "scopus", format = "bibtex")
s1 =  data.frame(convert2df("scopus_bottom2000F.bib", dbsource = "scopus", format = "bibtex"))
s2 =  data.frame(convert2df("scopus_top2000F.bib", dbsource = "scopus", format = "bibtex"))
common_cols <- intersect(colnames(s1), colnames(s2))
scopus  = rbind(
  subset(s1, select = common_cols), 
  subset(s2, select = common_cols)
)

scopus <- data.frame(scopus) %>% mutate_all(as.character())

scopus <- clean(scopus)

# # Supplement
# scopus <- scopus %>%
#   filter(DT == "ARTICLE")

scopus <- scopus %>%
  select(-AU_UN_NR)

cplot2 <- miss(
  scopus,
  scopus,
  fill = "#FFCD00FF", 
  filter = c("AB", "DI", "TI"),
  names = c("AB", "DI", "TI")
)

emf("scopus1.emf", width = 3, height = 2.5)
print(cplot2$p1)
dev.off()

emf("scopus2.emf", width = 3, height = 3.7)
print(cplot2$p2)
dev.off()

emf("scopus3.emf", width = 5.8, height = 4.5)
print(cplot2$p4)
dev.off()


emf("scopus4.emf", width = 5.8, height = 4.5)
print(cplot2$p4)
dev.off()

tiff("scopus4.tiff", width = 5.8 * 413, height = 4.5 * 413)
print(cplot2$p4)
dev.off()

## PubMed ###########################################################################

search <- 'COVID-19'

# pubmed_at_26062020 <- pmApiRequest(search, limit=27000, api_key = "") ## INSERT API key
# saveRDS(pubmed_at_26062020, "pubmed_at_26062020.RDS")

pubmed_at_26062020 <- readRDS("pubmed_at_26062020.RDS")
pm <- pmApi2df(pubmed_at_26062020, format = "bibliometrix")

pm <- pm %>% mutate_all(as.character())

pm <- clean(pm)

## Supplement
# pm <- pm %>%
#   filter(DT %in% c("JOURNAL ARTICLE", "CASE REPORTS"))

pm <- pm %>%
  select(-AU1_CO, -AU_CO)

cplot3 <- miss(
  pm,
  pm,
  fill = "#5C88DAFF", 
  filter = c("AB", "DI", "TI"),
  names = c("AB", "DI", "TI")
)

emf("pubmed1.emf", width = 3, height = 2.5)
print(cplot3$p1)
dev.off()

emf("pubmed2.emf", width = 3, height = 3.7)
print(cplot3$p2)
dev.off()

emf("pubmed3.emf", width = 5.8, height = 4.5)
print(cplot3$p4)
dev.off()

tiff("pubmed4.emf", width = 5.8 * 413, height = 2.9 * 413)
print(cplot3$p4)
dev.off()





### Compare
# Tocilizumab
topic1 <- "Tocilizumab"
# saveRDS(pubmedformating(topic1), paste(topic1, ".RDS", sep = ""))

pubmed1 <- readRDS(paste(topic1, ".RDS", sep = ""))
who1 <- read.csv(paste("w_", tolower(topic1), ".csv", sep = ""))
scopus1 <- read.csv(paste("s_", tolower(topic1), ".csv", sep = ""))

Venn(who1, pubmed1, scopus1, topic1)
data <- data.frame(Topic = topic1, who = nrow(who1), pubmed = nrow(pubmed1), scopus = nrow(scopus1))


l = matching(who1, scopus1, pubmed1, topic1)
topicx = topic1
xplot = plotting2(l)



saveplot(xplot, topicx)




# Depression

topic2 <- "Depression"
# saveRDS(pubmedformating(topic2), paste(topic2, ".RDS", sep = ""))

pubmed2 <- readRDS(paste(topic2, ".RDS", sep = ""))
who2 <- read.csv(paste("w_", tolower(topic2), ".csv", sep = ""))
scopus2 <- read.csv(paste("s_", tolower(topic2), ".csv", sep = ""))

Venn(who2, pubmed2, scopus2, topic2)
data <- updateData(who2, pubmed2, scopus2, topic2, data)

l = matching(who2, scopus2, pubmed2, topic2)
topicx = topic2
xplot = plotting2(l)
saveplot(xplot, topicx)


# Chloroquine

topic3 <- "Chloroquine"
# saveRDS(pubmedformating(topic3), paste(topic3, ".RDS", sep = ""))

pubmed3 <- readRDS(paste(topic3, ".RDS", sep = ""))
who3 <- read.csv(paste("w_", tolower(topic3), ".csv", sep = ""))
scopus3 <- read.csv(paste("s_", tolower(topic3), ".csv", sep = ""))

Venn(who3, pubmed3, scopus3, topic3)
data <- updateData(who3, pubmed3, scopus3, topic3, data)

l = matching(who3, scopus3, pubmed3, topic3)
topicx = topic3
xplot = plotting2(l)
saveplot(xplot, topicx)


# Azithromycin

topic4 <- "Azithromycin"
# saveRDS(pubmedformating(topic4), paste(topic4, ".RDS", sep = ""))

pubmed4 <- readRDS(paste(topic4, ".RDS", sep = ""))
who4 <- read.csv(paste("w_", tolower(topic4), ".csv", sep = ""))
scopus4 <- read.csv(paste("s_", tolower(topic4), ".csv", sep = ""))

Venn(who4, pubmed4, scopus4, topic4)
data <- updateData(who4, pubmed4, scopus4, topic4, data)

l = matching(who4, scopus4, pubmed4, topic4)
topicx = topic4
xplot = plotting2(l)
saveplot(xplot, topicx)


# Diabetes

topic5 <- "Diabetes"
# saveRDS(pubmedformating(topic5), paste(topic5, ".RDS", sep = ""))

pubmed5 <- readRDS(paste(topic5, ".RDS", sep = ""))
who5 <- read.csv(paste("w_", tolower(topic5), ".csv", sep = ""))
scopus5 <- read.csv(paste("s_", tolower(topic5), ".csv", sep = ""))

Venn(who5, pubmed5, scopus5, topic5)
data <- updateData(who5, pubmed5, scopus5, topic5, data)

l = matching(who5, scopus5, pubmed5, topic5)
topicx = topic5
xplot = plotting2(l)
saveplot(xplot, topicx)


# Hypertension

topic6 <- "Hypertension"
# saveRDS(pubmedformating(topic6), paste(topic6, ".RDS", sep = ""))

pubmed6 <- readRDS(paste(topic6, ".RDS", sep = ""))
who6 <- read.csv(paste("w_", tolower(topic6), ".csv", sep = ""))
scopus6 <- read.csv(paste("s_", tolower(topic6), ".csv", sep = ""))

Venn(who6, pubmed6, scopus6, topic6)
data <- updateData(who6, pubmed6, scopus6, topic6, data)

l = matching(who6, scopus6, pubmed6, topic6)
topicx = topic6
xplot = plotting2(l)
saveplot(xplot, topicx)


# Quarantine

topic7 <- "Quarantine"
# saveRDS(pubmedformating(topic7), paste(topic7, ".RDS", sep = ""))

pubmed7 <- readRDS(paste(topic7, ".RDS", sep = ""))
who7 <- read.csv(paste("w_", tolower(topic7), ".csv", sep = ""))
scopus7 <- read.csv(paste("s_", tolower(topic7), ".csv", sep = ""))

Venn(who7, pubmed7, scopus7, topic7)
data <- updateData(who7, pubmed7, scopus7, topic7, data)

l = matching(who7, scopus7, pubmed7, topic7)
topicx = topic7
xplot = plotting2(l)
saveplot(xplot, topicx)




# Shock
topic8 <- "Shock"
# saveRDS(pubmedformating(topic8), paste(topic8, ".RDS", sep = ""))

pubmed8 <- readRDS(paste(topic8, ".RDS", sep = ""))
who8 <- read.csv(paste("w_", tolower(topic8), ".csv", sep = ""))
scopus8 <- read.csv(paste("s_", tolower(topic8), ".csv", sep = ""))

Venn(who8, pubmed8, scopus8, topic8)
data <- updateData(who8, pubmed8, scopus8, topic8, data)

l = matching(who8, scopus8, pubmed8, topic8)
topicx = topic8
xplot = plotting2(l)
saveplot(xplot, topicx)


w = pubmed8
  w = who8
s = select(scopus, Abstract, Title, Index.Keywords, Author.Keywords, DOI)

for(i in 1:30){
  w[,i] = str_detect(toupper(w[,i]), toupper(topic))
}

w[is.na(w)]<-FALSE

w[w == TRUE]<-NA



pcompare <- tidyr::gather(data, cat, val, -Topic) %>% 
  ggplot(aes(x = Topic, y = val, fill = cat)) + 
    geom_col(position = "dodge") + 
    scale_fill_manual(values = c('#5C88DAFF', '#FFCD00FF',"#00B4E2FF")) +
    theme(legend.title = element_blank()) + 
    xlab("") + ylab("Number of articles [n]")


emf("compare.emf", width = 7.5, height = 3.4)
print(pcompare)
dev.off()




all_spec_who <- vis_miss(df %>% select(Abstract, Title, Descriptor.s.), cluster = TRUE) + theme(axis.text.x = element_text(angle = 0)) + coord_flip()
all_spec_sc <- vis_miss(scopus %>% select(AB, TI, DE, ID), cluster = TRUE) + theme(axis.text.x = element_text(angle = 0)) + coord_flip()
all_spec_pm <- vis_miss(pm %>% select(AB, TI, DE, ID), cluster = TRUE) + theme(axis.text.x = element_text(angle = 0)) + coord_flip()

emf(paste("st_", "all", "_pubmed.emf", sep = ""), width = 4.2, height = 2.5)
print(all_spec_pm)
dev.off()
emf(paste("st_", "all", "_scopus.emf", sep = ""), width = 4.2, height = 2.5)
print(all_spec_sc)
dev.off()
emf(paste("st_", "all", "_who.emf", sep = ""), width = 4.2, height = 2.5)
print(all_spec_who)
dev.off()


emf(paste("st_", "all", "_pubmedL.emf", sep = ""), width = 4.2*6, height = 2.5*6)
print(all_spec_pm)
dev.off()
emf(paste("st_", "all", "_scopusL.emf", sep = ""), width = 4.2*6, height = 2.5*6)
print(all_spec_sc)
dev.off()
emf(paste("st_", "all", "_whoL.emf", sep = ""), width = 4.2*6, height = 2.5*6)
print(all_spec_who)
dev.off()

tiff(paste("st_", "all", "_pubmed.tiff", sep = ""), width = 4.2*400, height = 2.5*400)
print(all_spec_pm)
dev.off()
tiff(paste("st_", "all", "_scopus.tiff", sep = ""), width = 4.2*400, height = 2.5*400)
print(all_spec_sc)
dev.off()
tiff(paste("st_", "all", "_who.tiff", sep = ""), width = 4.2*400, height = 2.5*400)
print(all_spec_who)
dev.off()







scopus_total <- read.csv("scopus_all.csv")

Venn(df, pm, scopus_total, "total")

## excluding abstracts

scopus_noab <- scopus[scopus$DT == "ARTICLE",] 
pm_noab <- pm[pm$DT %in% c("JOURNAL ARTICLE", "CASE REPORTS"),] 

table(is.na(scopus_noab$AB))
table(is.na(pm_noab$AB))






#####################
# WHO - 2020-05-19  #
#####################

df <- readxl::read_excel("CSV as at 19 May 2020-Full database.xlsx")


df <- clean(df)


cplot4 <- miss(
  df,
  df %>% mutate(Title.Authors.Journal = paste(Title, Authors, Journal)),
  fill = "#00B5E2FF", 
  filter = c("Abstract", "Accession Number", "DOI", "Title.Authors.Journal"),
  names = c("Abstract", "Accession\nNumber", "DOI", "Title\nAuthors\nJournal")
)


emf("whoold1.emf", width = 3, height = 2.5)
print(cplot4$p1)
dev.off()

emf("whoold2.emf", width = 3, height = 3.7)
print(cplot4$p2)
dev.off()

emf("whoold3.emf", width = 3, height = 2.7)
print(cplot4$p3)
dev.off()

emf("whoold4.emf", width = 5.8, height = 2.9)
print(cplot4$p4)
dev.off()

tiff("who4t.tiff", width = 5.8 * 413, height = 2.9 * 413)
print(cplot4$p4)
dev.off()

df$ID <- seq(1, dim(df)[1], 1)
df_doi <- df[is.na(df$DOI) != TRUE,] 
table(duplicated(df_doi$DOI) | duplicated(df_doi$DOI, fromLast=TRUE))
duplicateddoi <- df_doi[duplicated(df_doi$DOI) | duplicated(df_doi$DOI, fromLast=TRUE),]
final1 <- rbind(duplicateddoi)
final1 <- final1[duplicated(final1$ID) == FALSE,]
final1 %>%
  select(ID, Authors, Title, DOI, 'Accession Number', "Journal", Abstract) %>%
  arrange(DOI) %>%
  write.csv("WHO_older - duplicated based on DOI.csv", row.names = FALSE)


