# Coverage, inputs and overlap of the extended Unified Democracy Scores (democracyData::extended_uds)
# against the EIU baseline (di_score) in panel_data.csv. Run from the repo root.
# Result (2026-11-18): UDS runs 1502-2025; 159 countries x 8 yrs = 1272 matched country-years;
# EIU is an input in 100% of them; r = 0.95 pooled, 0.72 within-country, 0.50 on annual changes.

if (requireNamespace("democracyData", quietly=TRUE)) extended_uds <- democracyData::extended_uds else { tf <- tempfile(fileext=".rda"); download.file("https://raw.githubusercontent.com/xmarquez/democracyData/master/data/extended_uds.rda", tf, mode="wb"); load(tf) }; u <- as.data.frame(extended_uds[extended_uds$year>=2016 & extended_uds$year<=2023, c("extended_country_name","GWn","year","z1","num_measures","measures")])
p <- read.csv("data/output_datasets/panel_data.csv")
pc <- unique(p[,c("country_name","country_code")])
norm <- function(x) gsub("[^a-z]","",iconv(tolower(gsub("\\s*\\(.*\\)","",x)),to="ASCII//TRANSLIT"))
un <- unique(u$extended_country_name)
pc$m <- un[match(norm(pc$country_name), norm(un))]
man <- c(BHS="Bahamas",BLZ="Belize",BRN="Brunei",BIH="Bosnia-Herzegovina",CIV="Cote D'Ivoire",COD="Congo, Democratic Republic of (Zaire)",COG="Congo",
 DEU="German Federal Republic",IRN="Iran (Persia)",ITA="Italy/Sardinia",KOR="Korea, Republic of",MKD="Macedonia (Former Yugoslav Republic of)",
 STP="Sao Tome and Principe",SUR="Surinam",SWZ="Swaziland",TTO="Trinidad and Tobago",TZA="Tanzania/Tanganyika",USA="United States of America",
 VNM="Vietnam, Democratic Republic of",YEM="Yemen (Arab Republic of Yemen)",CZE="Czech Republic")
for(k in names(man)) pc$m[pc$country_code==k] <- man[k]
cat("still unmatched:\n"); print(pc[is.na(pc$m),1:2])
# check dup matches
print(pc[duplicated(pc$m)&!is.na(pc$m),])
p$m <- pc$m[match(p$country_code,pc$country_code)]
d <- merge(p[,c("country_code","m","year","di_score")], u, by.x=c("m","year"), by.y=c("extended_country_name","year"), all.x=TRUE)
cat("panel rows:",nrow(p)," with di_score:",sum(!is.na(p$di_score))," with uds:",sum(!is.na(d$z1))," both:",sum(!is.na(d$z1)&!is.na(d$di_score)),"\n")
cat("countries with both:", length(unique(d$country_code[!is.na(d$z1)&!is.na(d$di_score)])),"\n")
b <- d[!is.na(d$z1)&!is.na(d$di_score),]
cat("pooled r =", round(cor(b$di_score,b$z1),3), " spearman =", round(cor(b$di_score,b$z1,method="spearman"),3),"\n")
for(y in 2016:2023){x<-b[b$year==y,]; cat(y,nrow(x),round(cor(x$di_score,x$z1),3),"\n")}
# within-country (demeaned) correlation
dm <- function(v,g) v-ave(v,g)
cat("within-country r =", round(cor(dm(b$di_score,b$country_code),dm(b$z1,b$country_code)),3),"\n")
# year-to-year changes
b <- b[order(b$country_code,b$year),]
b$ddi <- ave(b$di_score,b$country_code,FUN=function(v)c(NA,diff(v))); b$duds <- ave(b$z1,b$country_code,FUN=function(v)c(NA,diff(v)))
cat("r of annual changes =", round(cor(b$ddi,b$duds,use="complete"),3),"\n")
cat("share of panel country-years where EIU is a UDS input:", round(mean(sapply(b$measures,function(m)"eiu"%in%m)),3),"\n")
cat("mean num measures in panel country-years:", round(mean(b$num_measures),2),"; by year:\n"); print(tapply(b$num_measures,b$year,mean))
cat("measures seen in panel 2016-2023:\n"); print(sort(table(unlist(b$measures)),decreasing=TRUE))
cat("measures 2023:\n"); print(sort(table(unlist(b$measures[b$year==2023])),decreasing=TRUE))
