#Open City_Country and PM10 PM2.5 tables
Readings <- read.tidy("WHO_data_average_2013_2016_subset.txt")
Country <- read.tidy("WHO_data.txt")[c("city","country")]
Country <- unique(Country)

Merge <- merge(Readings,Country,by.x="City", by.y="city",all.x=TRUE)

Merge <- unique(Merge)

Merge$country <- factor(Merge$country, levels = c("Australia","Brazil","Canada","China","France","Japan","United Kingdom","United States of America","Baoding","Dalian"))

#Plot
Plot <- ggplot(Merge,aes(x=PM10,y=PM25, color = country))
Plot <- Plot + geom_point(size=2) + theme_classic() + scale_color_brewer(palette="Set3")
ggsave("WHO_average_2013_2016.pdf",width=7,height=9,units="in")
