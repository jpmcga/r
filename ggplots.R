library(ggplot2)



gg <- ggplot(diamonds, aes(x=carat, y=price, color=cut)) +
  geom_point() +
  labs(title = "W.e", x="Carat", y="Price")

gg1 <- gg +
  theme(plot.title=element_blank(), 
        axis.text.x=element_text(size=10), 
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=25),
        axis.title.y=element_text(size=25),
        legend.position = "none") + 
  facet_grid(color ~ cut)
gg1

