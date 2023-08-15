library(readxl)
library(tidyverse)

AL_est<-read_excel('~/Downloads/41586_2018_836_MOESM3_ESM.xlsx')

params_inferred <- tibble(r_on = numeric(),
                          r_off = numeric(),
                          r_prod = numeric())
for (i in 1:nrow(AL_est)) {
  
  tmp <- substr(as.character(AL_est[i,2]), 2,nchar(as.character(AL_est[i,2]))-2)
  tmp2 <- as.numeric(strsplit(tmp, split = ' ')[[1]])
  tmp3 <- tmp2[!is.na(tmp2)]
  tmp4 <- tibble(
    r_on = tmp3[1],
    r_off = tmp3[2],
    r_prod = tmp3[3]
  )
  
  params_inferred %<>% bind_rows(tmp4)
  
}

params_inferred %<>%
  mutate(on_off_ratio = r_on/r_off)

p1<-ggplot(params_inferred, aes(x=r_on)) + geom_histogram() + theme_classic()
p2<-ggplot(params_inferred, aes(x=log10(on_off_ratio))) + geom_histogram() + theme_classic()
p3<-ggplot(params_inferred, aes(x=r_prod)) + geom_histogram() +scale_x_log10() + theme_classic()

svglite('~/code/grn_nitc/larsson/parsed_params.svg', width = 6, height = 4)
grid.arrange(p1,p2,p3, ncol = 3)
dev.off();

ggplot() + geom_histogram(data = params_inferred %>% pivot_longer(cols = r_on:on_off_ratio, names_to = 'parameter', values_to = 'value'), aes(value)) + facet_wrap(~parameter, scales = 'free_x')
