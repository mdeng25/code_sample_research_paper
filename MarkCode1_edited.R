mRNAvRiboDiff <- dat.gather.rpm.differences %>% dplyr::select(Geneid, WT_D_m, WT_D)
mRNAvRiboDiff.distance <- mutate(mRNAvRiboDiff, distanceFromxy = (abs(WT_D_m - WT_D))/(sqrt(2)))
# mRNAvRiboDiff.distance <- merge(mRNAvRiboDiff.distance, dat.gather.rpm, by.x = 'Geneid')
mRNAvRiboDiff.distance2 <- tbl_df(merge(mRNAvRiboDiff.distance, dat.gather.rpm, by.x = 'Geneid'))
mRNAvRiboDiff.distance3 <- mRNAvRiboDiff.distance2 %>% dplyr::select(Geneid, ALIAS, WT_D_m, WT_D, distanceFromxy, Length, rpm, sample)
mRNAvRiboDiff.distance4 <- mRNAvRiboDiff.distance3[(mRNAvRiboDiff.distance3$sample=="WT_D"),]
mRNAvRiboDiff.distance5 <- mutate(mRNAvRiboDiff.distance4, RPKM = rpm/(Length/1000))

#plot discrepancy vs gene length no log
plot.discrepancyvslength <- ggplot(mRNAvRiboDiff.distance5, aes(x=Length, y=distanceFromxy)) + geom_point() + scale_x_continuous(limits=c(0,1000))

## Name the above plot and save it

## Instead of taking Log_10 of Length, try restricting the limits of the x-axis such that outliers don't dominate
## Add: + scale_x_continuous(limits=c(0,1000))

#plot discrepency vs gene length, log 10 on x-axis
#plot.discrepancyvslength.log10 <- ggplot(mRNAvRiboDiff.distance5, aes(x=log(Length,10), y=distanceFromxy)) + geom_point() #+ scale_y_continuous(limits=c(0,1000))

## Look in GeneSets.R for a well-annotated gene sets

#plot ribosome and proteasome
mRNAvRiboDiff.distance.rp <- mRNAvRiboDiff.distance5 %>% filter(ALIAS %in% ribo.symbols) 
mRNAvRiboDiff.distance.psm <- mRNAvRiboDiff.distance5 %>% filter(ALIAS %in% proteasome)

plot.rp <- ggplot(mRNAvRiboDiff.distance.rp, aes(x=ALIAS, y=RPKM)) + geom_point() +  geom_text(aes(label=ALIAS),hjust=0, vjust=0)

plot.psm <- ggplot(mRNAvRiboDiff.distance.psm, aes(x=ALIAS, y=RPKM)) + geom_point() + geom_text(aes(label=ALIAS),hjust=0, vjust=0)
