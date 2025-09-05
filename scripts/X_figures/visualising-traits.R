input <- '.'
library(ggplot2)
# --> LOAD TRAITS 
Tr <- read.csv(file.path(input,"data/1_preprocessing/Tr_aits/traits-guild_migration.csv"),row.names = 2)[,c(2,3)]
Tr <- Tr[sort(row.names(Tr)),]

Tr$Migration_a3_DOF <- factor(Tr$Migration_a3_DOF,
                              level = c('sedentary',
                                        'sedentary and short-distance',
                                        'short-distance',
                                        'short-and long-distance',
                                        'long-distance'))

pdf(file.path(input,'figs/1_preprocessing/guilds-strategies/guilds-strategies.pdf'),
    width = 10,height=6)
mig_freq <- table(Tr$Migration_a3_DOF)

ggplot(Tr,aes(Migration_a3_DOF))+
  geom_bar(na.rm = T) +
  labs(title = paste0('Nr of species per migratory strategy (n = ',nrow(Tr),')'),
       y = 'Nr of species',
       x = NULL) +
  scale_x_discrete(guide = guide_axis(n.dodge=2))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14),
        plot.margin = unit(c(1,1,1,1),'cm'))


guild_freq <- table(Tr$foraging_guild_consensus)
order(guild_freq)


ggplot(Tr,aes(foraging_guild_consensus))+
  geom_bar(na.rm = T) +
  labs(title = paste0('Nr of species per foraging guild (n = ',nrow(Tr),')'),
       y = 'Nr of species',
       x = NULL) +
  scale_x_discrete(guide = guide_axis(angle=45))+
  theme(axis.text = element_text(size=11),
        axis.title = element_text(size=14),
        plot.margin = unit(c(1,1,1,1),'cm'))

dev.off()
