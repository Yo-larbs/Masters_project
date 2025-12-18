## ----setup, include=FALSE--------------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  dev = ifelse(knitr::is_latex_output(), "pdf", "svg"),
  fig.ext = ifelse(knitr::is_latex_output(), "pdf", "svg"),
  fig.path = "../../output/figures/cleaning_and_classification/virus/"
  # write files with .svg extension
  )


## --------------------------------------------------------------------------------------------------------
#install.packages("factoextra")
#install.packages("tidyr")
#install.packages("readr")
#install.packages("ggplot2")
#install.packages("cellranger")
#install.packages('ggrepel')
#install.packages("minpack.lm")
#install.packages("drc")
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(ggrepel)
library(minpack.lm)
library(drc)
library(factoextra)


## --------------------------------------------------------------------------------------------------------
setwd('~/Documents/Masters_project/')
HIV1_virus_data=read_csv('inputs/Cleaning_and_classification/hiv1_virus_collated_n0_uncleaned(Sheet1).csv')


## --------------------------------------------------------------------------------------------------------
HIV1_virus_data=HIV1_virus_data[,c(2,6:9,11:15,17)]
HIV1_virus_data=na.omit(HIV1_virus_data)
HIV1_virus_data=HIV1_virus_data%>%mutate(cell_line=if_else(cell_line=="Blix_1","Bilx_1",HIV1_virus_data$cell_line))



## --------------------------------------------------------------------------------------------------------
library(dplyr)
HIV1_virus_data=HIV1_virus_data%>%group_by(cell_line,screen_nb,replicate)%>% 
  mutate(max_gfp=max(assay_output),max_gfp_titre=which.max(assay_output)) %>%
mutate(assay_output=if_else((infection_volume_ul>infection_volume_ul[max_gfp_titre] & assay_output  <= 0.8 * max_gfp),max_gfp,assay_output))%>%
  dplyr::select(-max_gfp,-max_gfp_titre) %>% 
  ungroup()


## --------------------------------------------------------------------------------------------------------
virus_data_per_cell_line=HIV1_virus_data%>%nest_by(cell_line,.keep = T)


## --------------------------------------------------------------------------------------------------------
apply_plot=function(i){ggplot(data = i,aes(x=infection_volume_ul,y=assay_output,group = interaction(screen,replicate) ))+
  geom_point()+
  ggtitle(i$cell_line[1])+
  theme_classic()
}


## ----scatter plot----------------------------------------------------------------------------------------
scatter_plot_virus=purrr::map(virus_data_per_cell_line$data,apply_plot)
names(scatter_plot_virus)=virus_data_per_cell_line$cell_line
scatter_plot_virus


## --------------------------------------------------------------------------------------------------------
apply_plot_boxplot=function(i){ggplot(data = i,aes(y=assay_output,x=titre,))+
  geom_boxplot()+
  ggtitle(i$cell_line[1])+
  #stat_summary(geom="text", fun.y=quantile,
               #aes(label=sprintf("%1.1f", ..y..), color=factor(titre)),
               #position=position_nudge(x=0.6), size=3.5) +
  theme_classic()
}


## ----boxplot---------------------------------------------------------------------------------------------
box_plot_virus=purrr::map(virus_data_per_cell_line$data,apply_plot_boxplot)
names(box_plot_virus)=virus_data_per_cell_line$cell_line
box_plot_virus


## --------------------------------------------------------------------------------------------------------
class(HIV1_virus_data)


## ----unnormalised_boxplot--------------------------------------------------------------------------------
unnormalised_boxplot=ggplot(data = HIV1_virus_data,aes(y=assay_output, x=as.factor(screen_nb),group = screen_nb,))+
  geom_boxplot()+
  ggtitle('Before normalisation')+
  xlab("Screen")+
  ylab("Supernatant Gag")+
  theme_classic()
bmp("unnormalised_boxplot.bmp",width = 3200 ,height = 2400,res = 500,pointsize = 300)
unnormalised_boxplot
dev.off()


## --------------------------------------------------------------------------------------------------------

HIV1_virus_data=HIV1_virus_data%>%
  group_by(screen_nb)%>%
  mutate(zscore=((((assay_output-mean(assay_output))/sd(assay_output)))))%>%
  ungroup()

virus_data_per_cell_line=HIV1_virus_data%>%nest_by(cell_line,.keep = T)



## ----normalised_boxplot----------------------------------------------------------------------------------
normalised_boxplot=ggplot(data = HIV1_virus_data,aes(y=zscore, x=as.factor(screen_nb),group = screen_nb,))+
  geom_boxplot()+
  ggtitle('after normalisation')+
  xlab("Screen")+
  ylab("Zscore")+
  theme_classic()
bmp("normalised_boxplot.bmp",width = 3200 ,height = 2400,res = 500,pointsize = 300)
normalised_boxplot
dev.off()


## --------------------------------------------------------------------------------------------------------
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}



## --------------------------------------------------------------------------------------------------------
outlier_function=function(i){
  i %>%
  group_by(titre) %>%
  mutate(outlier = ifelse(is_outlier(zscore), zscore, as.numeric(NA))) 
}



## --------------------------------------------------------------------------------------------------------
all_outliers=purrr::map(virus_data_per_cell_line$data,outlier_function)

HIV1_virus_data_outliers=bind_rows(all_outliers)%>%
  group_by(cell_line,screen_nb,replicate)%>%
  filter(((sum(is.na(outlier)))<3))

HIV1_virus_data=bind_rows(all_outliers)%>%group_by(cell_line,screen_nb,replicate)%>%filter(((sum(is.na(outlier)))>2))



## --------------------------------------------------------------------------------------------------------
cell_for_plotting=c("CTR_M2_O5","TZM","Vuud_2")

df=(bind_rows(all_outliers))[(bind_rows(all_outliers)$cell_line)%in%unique(HIV1_virus_data_outliers$cell_line[!HIV1_virus_data_outliers$cell_line%in%cell_for_plotting]),]%>%group_by(cell_line,screen_nb,replicate)%>%mutate(max_assay=max(zscore))

df2=HIV1_virus_data_outliers[!HIV1_virus_data_outliers$cell_line%in%cell_for_plotting,]%>%group_by(cell_line,screen_nb,replicate)%>%mutate(max_assay=max(zscore))


ggplot()+
  geom_line(data=df,aes(y=(max_assay),x = interaction(replicate,screen),group = 1))+
  geom_point(data=df2,aes(y=(max_assay),x = interaction(replicate,screen)),colour="red")+
    geom_point(data=df,aes(y=(zscore),x = interaction(replicate,screen)),colour="black",shape=2,size=0.2)+

  facet_wrap(~cell_line, scales = "free_y") +
  #geom_boxplot(data=df,aes(y=(zscore),x = interaction(replicate,screen),),show.legend = T)+
  
  theme_classic()+
  #ggtitle()+
  xlab("Replicate.Screen")+
  ylab("Ζscore")+
  theme(text = element_text(size = 10))



## --------------------------------------------------------------------------------------------------------
HIV1_virus_data=HIV1_virus_data%>%group_by(cell_line,screen,screen_nb,replicate)%>%
  mutate(area_under_curve=AUC(infection_volume_ul,(zscore)))%>%
  ungroup()


## --------------------------------------------------------------------------------------------------------
virus_data_per_cell_line=HIV1_virus_data%>%nest_by(cell_line,.keep = T)


## --------------------------------------------------------------------------------------------------------
logarithmic_func <- function(x, a, b) {
  a + b * log((x),base = 4)
}



## --------------------------------------------------------------------------------------------------------
fit_log <- function(data) {

  # Fit the logarithmic model for each group and add coefficients to the data frame
  fitted_data <- data %>%
    group_by(screen_nb,replicate) %>%
    group_modify(~ {
      tryCatch({
        # Fit the logarithmic model
        model <- nlsLM(
          zscore ~ logarithmic_func(infection_volume_ul, a, b),#, c),
          data = .x,
          start = list(a = 0, b = 1),#, c = -0.0001),
          control = nls.control(maxiter = 1024)
        )
        
        # Extract coefficients
        coefs <- coef(model)
        
        # Add coefficients back to the group's data
        .x %>%
          mutate(a = coefs["a"], b = coefs["b"])#, c = coefs["c"])
      }, error = function(e) {
        # If fitting fails, add NA values for coefficients
        .x %>%
          mutate(a = NA, b = NA)#, c = NA)
      })
    }) %>%
    ungroup()  # Ungroup after processing
  
  return(fitted_data)
}
# Assuming your data frame is named 'df' and has columns 'screen', 'infection_volume_ul', and 'mean(zscore)'
virus_data_per_cell_line=purrr::map(virus_data_per_cell_line$data,fit_log)%>%bind_rows()%>%nest_by(cell_line,.keep = T)



## --------------------------------------------------------------------------------------------------------
apply_plot_log=function(i){
i <- i %>%
    mutate(predicted = logarithmic_func(infection_volume_ul, i$a, i$b))#,i$c))
  
  # Create the plot
  p=ggplot(data = i, aes(x = infection_volume_ul, y = zscore, group = interaction(screen_nb,replicate))) +
    geom_point(aes(color=screen_nb),color = "blue", size = 2, alpha = 0.6) + # Actual data points
    geom_line(aes(y = predicted,color=screen_nb), color = "red", size = 1) + # Fitted line
    ggtitle(unique(i$cell_line)) + # Dynamic title
    labs(x = "Infection Volume (uL)", y = "Mean Z-Score") +
    theme_classic()
  
    params_text <- paste0(
    "Parameters:\n",
    "a = ", round(mean(i$a, na.rm = TRUE), 2), "\n",
    "b = ", round(mean(i$b, na.rm = TRUE), 2), "\n"
    #"d = ", round(mean(i$c, na.rm = TRUE), 2), "\n"
  )
  
  p + annotate("text",
               x = Inf, y = -Inf,
               label = params_text,
               hjust = 1.1, vjust = -0.1,
               size = 3)
}


## ----logarithmic-----------------------------------------------------------------------------------------
logarithmic_plot_virus=purrr::map(virus_data_per_cell_line$data,apply_plot_log)
names(logarithmic_plot_virus)=virus_data_per_cell_line$cell_line
logarithmic_plot_virus


## --------------------------------------------------------------------------------------------------------
# Correct 4PL formula (parentheses fix)
logistic_func <- function(x, b, d, e,c){
 c+ ((d-c)/(1 + exp(b*(x - e))))# Fixed denominator placement
}


logistic_fit <- function(data) {
    fitted_data <- data %>%
        group_by(screen_nb,replicate) %>% 
        group_modify(~ {
            tryCatch({ 
                model <- drm(zscore ~ infection_volume_ul, fct = L.4(), data = .x)
                coefs <- coef(model)

                .x %>%
                    mutate(
                        logis_b = coefs["b:(Intercept)"],
                        logis_e = coefs["e:(Intercept)"],
                        logis_c = coefs["c:(Intercept)"],
                        logis_d = coefs["d:(Intercept)"]
                    )
            }, error = function(e) {
                warning(paste("drc error for screen_nb:", .x$screen_nb[1], ":", e$message))
                .x %>%
                    mutate(logis_b = NA, logis_e = NA, logis_d = NA)
            })
        }) %>%
        ungroup()
    return(fitted_data)
}

# Apply to cell line data
virus_data_per_cell_line <- purrr::map(virus_data_per_cell_line$data,logistic_fit)%>%bind_rows()%>%nest_by(cell_line,.keep = T)




## --------------------------------------------------------------------------------------------------------
apply_plot_logistic <- function(i) {
  # Generate predictions using row-wise parameters
  plot_data <- i %>%
    mutate(
      predicted = logistic_func(
        x = infection_volume_ul,
        b = logis_b,
        c = logis_c,
        d = logis_d,
        e = logis_e
      )
    )
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = infection_volume_ul, y = zscore)) +
    # Actual data points
    geom_point(
      aes(color = factor(screen_nb)),  # Color by screen_nb
      size = 2,
      alpha = 0.6,
      show.legend = T
    ) +
    # Fitted curves
    geom_line(
      aes(y = predicted, group = interaction(screen_nb, replicate)),
      color = "red",
      linewidth = 1
    ) +
    # Labels and theme
    labs(
      title = unique(plot_data$cell_line),
      x = "Infection Volume (μL)",
      y = "Mean Z-Score"
    ) +
    theme_classic(base_size = 12)
    #scale_x_log10()  # Often useful for dose-response data
  
  # Add regression parameters annotation
  params_text <- paste0(
    "Parameters:\n",
    "b = ", round(mean(plot_data$logis_b, na.rm = TRUE), 2), "\n",
    "c = ", round(mean(plot_data$logis_c, na.rm = TRUE), 2), "\n",
    "d = ", round(mean(plot_data$logis_d, na.rm = TRUE), 2), "\n",
    "e = ", round(mean(plot_data$logis_e, na.rm = TRUE), 2)
  )
  
  p + annotate("text",
               x = Inf, y = -Inf,
               label = params_text,
               hjust = 1.1, vjust = -0.1,
               size = 3)
}


## ----logistic--------------------------------------------------------------------------------------------
logistic_plot_virus=purrr::map(virus_data_per_cell_line$data,apply_plot_logistic)
names(logistic_plot_virus)=virus_data_per_cell_line$cell_line
logistic_plot_virus


## --------------------------------------------------------------------------------------------------------
SD_table=data.frame()

for (i in virus_data_per_cell_line[2][[1]]){
  SD_table=rbind(SD_table,distinct(i%>%group_by(screen_nb,replicate)%>% 
                                     mutate(max_zscore=max(zscore))%>%
                                     ungroup()%>%
                                     reframe(
                                       cell_line=cell_line,
                                       mean_max_zscore=mean(max_zscore),
                                       SD_max_zscore=sd(max_zscore),
                                       CV_z=(SD_max_zscore/mean_max_zscore)*100,
                                       mean_AUC=mean(abs(area_under_curve)),
                                       SD_AUC=sd((abs(area_under_curve))),
                                       CV_AUC  = (SD_AUC/mean_AUC)*100)
                                   )
                 )
}


## --------------------------------------------------------------------------------------------------------
mean(SD_table$SD_max_zscore[c(1:10,12:113,115:152)])
mean(na.omit(SD_table$SD_AUC))#[c(1:10,12:113,115:152)])


## --------------------------------------------------------------------------------------------------------
HIV1_virus_data=bind_rows(virus_data_per_cell_line)%>%unnest()

HIV1_virus_data_mean=HIV1_virus_data%>%group_by(batch,cell_line,screen,screen_nb,titre,infection_volume_ul)%>%
  reframe(assay_output=mean(assay_output),
            zscore=mean(zscore),
            a_log=mean(a),
            b_log=mean(b),
            logis_b=mean(logis_b),
            logis_d=mean(logis_d),
            logis_c=mean(logis_c),
            logis_e=mean(logis_e),
            area_under_curve=mean(area_under_curve))

HIV1_virus_data_mean <- HIV1_virus_data_mean %>%
  mutate(across(c(zscore,
                  a_log, 
                  b_log,
                  logis_b, 
                  logis_d,
                  logis_c,
                  logis_e,
                  area_under_curve), as.numeric))



## --------------------------------------------------------------------------------------------------------
HIV1_virus_data_PCA=unique(HIV1_virus_data_mean%>%group_by(cell_line,screen_nb)%>%
                             reframe(assay_output=max(zscore),
                                       a_log=a_log,
                                       b_log=(b_log),
                                       logis_b=logis_b,
                                       logis_d=logis_d,
                                       logis_c=logis_c,
                                       logis_e=logis_e,
                                       area_under_curve=area_under_curve))


## --------------------------------------------------------------------------------------------------------

HIV1_virus_data_PCA_sum <- HIV1_virus_data_PCA %>%
  group_by(cell_line) %>%
  reframe(
    assay_output      = mean(assay_output, na.rm = TRUE),
    a_log             = mean(a_log, na.rm = TRUE),
    b_log             = mean(b_log, na.rm = TRUE),
    #logis_b           = mean(logis_b, na.rm = TRUE),
    #logis_d           = mean(logis_d, na.rm = TRUE),
    #logis_c           = mean(logis_c, na.rm = TRUE),
    #logis_e           = mean(logis_e, na.rm = TRUE),
    area_under_curve  = abs(mean(area_under_curve, na.rm = TRUE))
  )


## --------------------------------------------------------------------------------------------------------
HIV1_virus_data_PCA_1=HIV1_virus_data_PCA_sum%>%mutate_if(is.numeric,scale)


HIV1_virus_data_PCA_1=subset(HIV1_virus_data_PCA_1, select=-cell_line)
HIV1_virus_data_PCA_1=as.matrix(HIV1_virus_data_PCA_1)
row.names(HIV1_virus_data_PCA_1)=HIV1_virus_data_PCA_sum$cell_line





## --------------------------------------------------------------------------------------------------------
HIV1_virus_data_PCA_1=as.data.frame(HIV1_virus_data_PCA_1)
dim(HIV1_virus_data_PCA_1)
#cor(x=HIV1_virus_data_PCA_1$c_log,y=HIV1_virus_data_PCA_1$area_under_curve)


## ----wssplot---------------------------------------------------------------------------------------------
wssplot <- function(data, nc=15, seed=1234){
                  wss <- (nrow(data)-1)*sum(apply(data,2,var))
                      for (i in 2:nc){
                set.seed(seed)
                    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
              plot(1:nc, wss, type="b", xlab="Number of Clusters",
                            ylab="Within groups sum of squares")
              wss
}
wssplot(HIV1_virus_data_PCA_1)


## --------------------------------------------------------------------------------------------------------
PCA=prcomp(x = HIV1_virus_data_PCA_1)


## ----scree_plot------------------------------------------------------------------------------------------
library(factoextra)
fviz_eig(PCA)


## ----PCA_spread------------------------------------------------------------------------------------------
PCA_spread=fviz_pca_ind(PCA,col.ind = "cos2",repel = F)
PCA_spread


## ----PCA contributions-----------------------------------------------------------------------------------
contributions_PCA=fviz_pca_var(PCA,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE ) 
#bmp("contributions_PCA_virus.bmp",width = 3200,height = 2400,res = 500)
contributions_PCA
#dev.off()


## ----biPlot----------------------------------------------------------------------------------------------
PCA_biplot=fviz_pca_biplot(PCA, repel=F,
                col.var = "#2E9FDF", # Variables color
                col.ind = "#696969"  # Individuals color
                )
PCA_biplot


## ----PCA dim 1-------------------------------------------------------------------------------------------
res_var_virus=get_pca_var(PCA)
res.ind <- get_pca_ind(PCA)
dim(res.ind$coord)
ind_coord_virus <- as.data.frame(res.ind$coord)
ind_coord_virus=ind_coord_virus[-which(row.names(ind_coord_virus)==('TZM')|row.names(ind_coord_virus)=='CTR_M2_O5'|row.names(ind_coord_virus)=='CTR_M2_O5_21'),]
res.ind$cos2=res.ind$cos2[-which(row.names(res.ind$cos2)==('TZM')|row.names(res.ind$cos2)=='CTR_M2_O5'|row.names(res.ind$cos2)=='CTR_M2_O5_21'),]


PCA_dim1=ggplot(data = ind_coord_virus,aes(x=Dim.1,y = 0,label=row.names(ind_coord_virus)))+
  geom_point(size = 3,aes(colour = res.ind$cos2[,1]), position = position_jitter(height = 0.02,width = 0),)+
  geom_label_repel(aes(label=ifelse(Dim.1 < quantile(Dim.1, 0.1) | Dim.1 > quantile(Dim.1, 0.9),row.names(ind_coord_virus),'')),max.overlaps = 150)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray")+
  theme_minimal() +
  xlab("PC1")+
  coord_cartesian(ylim = c(-0.05, 0.05))+
  theme(axis.ticks.y = element_blank(),
        axis.text.y  = element_blank(),
        panel.grid.major.y = element_blank())
         
PCA_dim1


## ----PCA classification----------------------------------------------------------------------------------
susceptible_virus_cell=ind_coord_virus[(ind_coord_virus$Dim.1 < quantile(ind_coord_virus$Dim.1, 0.1)),]
Resistant_virus_cell=ind_coord_virus[(ind_coord_virus$Dim.1 > quantile(ind_coord_virus$Dim.1, 0.9)),]

HIV_virus_data_susceptible=HIV1_virus_data_PCA_1[match(row.names(susceptible_virus_cell),(HIV1_virus_data_PCA_sum$cell_line)),]
HIV_virus_data_resistant=HIV1_virus_data_PCA_1[match(row.names(Resistant_virus_cell),(HIV1_virus_data_PCA_sum$cell_line)),]


set.seed(39)
clusters <- kmeans(ind_coord_virus[, c("Dim.1", "Dim.2")], centers = 10)


ind_coord_virus$group <- as.factor(clusters$cluster)
ind_coord_virus$pc_combined=0.80*ind_coord_virus$Dim.1 + 0.20*ind_coord_virus$Dim.2

PCA_classification=ggplot(data = ind_coord_virus,aes(x=Dim.1,y = Dim.2,label=row.names(ind_coord_virus)))+
  geom_point(size = 3,aes(colour = group))+
  geom_label_repel(aes(label=ifelse(Dim.1 < quantile(Dim.1, 0.1) | Dim.1 > quantile(Dim.1, 0.9),row.names(ind_coord_virus),'')),max.overlaps = 150)+
  xlab("PC1")+
  ylab("PC2")+
  theme_classic()


susceptible_virus_cell=ind_coord_virus[(ind_coord_virus$pc_combined < quantile(ind_coord_virus$pc_combined, 0.1)),]
Resistant_virus_cell=ind_coord_virus[(ind_coord_virus$pc_combined > quantile(ind_coord_virus$pc_combined, 0.9)),]



#bmp("PCA_classification.bmp",width = 3200,height = 2400,res = 500)
PCA_classification
#dev.off()

HIV1_virus_data_notable <- as.data.frame(rbind(HIV_virus_data_resistant, HIV_virus_data_susceptible)) %>%
  mutate(phenotype = ifelse(row.names(.) %in% row.names(HIV_virus_data_resistant),
                            "resistant",
                            "susceptible"))


## --------------------------------------------------------------------------------------------------------
susceptible_virus_clusters=ind_coord_virus[(ind_coord_virus$group==9)|(ind_coord_virus$group==5),]
resistant_virus_clusters=ind_coord_virus[(ind_coord_virus$group==8),]#|(ind_coord_virus$group==6),]

HIV_virus_data_susceptible=HIV1_virus_data_PCA_1[match(row.names(susceptible_virus_clusters),(HIV1_virus_data_PCA_sum$cell_line)),]
HIV_virus_data_resistant=HIV1_virus_data_PCA_1[match(row.names(resistant_virus_clusters),(HIV1_virus_data_PCA_sum$cell_line)),]

HIV1_virus_data_notable <- as.data.frame(rbind(HIV_virus_data_resistant, HIV_virus_data_susceptible)) %>%
  mutate(phenotype = ifelse(row.names(.) %in% row.names(HIV_virus_data_resistant),
                            "resistant",
                            "susceptible"))


## ----Z score split---------------------------------------------------------------------------------------
Z_score_split=ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=assay_output,colour = phenotype))+
  geom_point()+
   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
Z_score_split


## ----a_log split-----------------------------------------------------------------------------------------
a_log_split=ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=a_log,colour = phenotype))+
  geom_point()+
   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
a_log_split


## ----b_log split-----------------------------------------------------------------------------------------
b_log_split=ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=b_log,colour = phenotype))+
  geom_point()+
   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
b_log_split


## ----AUC_split-------------------------------------------------------------------------------------------
AUC_split=ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=area_under_curve,colour = phenotype))+
  geom_point()+
   theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
AUC_split


## ----b_logis---------------------------------------------------------------------------------------------
#b_logis=ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=logis_b,colour = phenotype))+
  #geom_point()


## --------------------------------------------------------------------------------------------------------
#ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=logis_d,colour = phenotype))+
  #geom_point()


## --------------------------------------------------------------------------------------------------------
#ggplot(HIV1_virus_data_notable,aes(x=row.names(HIV1_virus_data_notable),y=logis_e,colour = phenotype))+
  #geom_point()


## --------------------------------------------------------------------------------------------------------
virus_ranking=list()
for(col in names(HIV1_virus_data_PCA_1)) {
  
  # Extract the column values
  vec <- HIV1_virus_data_PCA_1[[col]]
  
  # Calculate the 10th and 90th quantiles (ignoring NAs)
  low_quant  <- quantile(vec, 0.1, na.rm = TRUE)
  high_quant <- quantile(vec, 0.9, na.rm = TRUE)
  
  # Identify rows where the value is either below the 10th percentile or above the 90th percentile
  extreme_idx <- which(vec < low_quant | vec > high_quant)
  
  # Subset the original data frame for these extreme values
  subset_data <- HIV1_virus_data_PCA_1[extreme_idx, ]
  
  # Add a phenotype column: if the value is below the 10th percentile, call it "resistant",
  # if above the 90th percentile, label it "susceptible"
  # (This uses the current column's values from the subset.)
  subset_data$phenotype <- ifelse(vec[extreme_idx] < low_quant, "resistant", "susceptible")
  subset_data$cell_line=rownames(subset_data)
  rownames(subset_data) <- NULL
  
  # Store the subset in the virus_ranking list using the column name as the key
  virus_ranking[[col]] <- subset_data
}





## --------------------------------------------------------------------------------------------------------
virus_ranking_df=bind_rows(virus_ranking)
cell_line_counts <- virus_ranking_df %>% 
  count(cell_line,phenotype)

# 3. Get the cell lines that appear more than 9 times.
common_virus_cell_lines <- cell_line_counts[cell_line_counts$n>3,]

# 4. Filter the original data to include only those common_virus cell lines.
most_common_virus <- inner_join(virus_ranking_df,common_virus_cell_lines,join_by(x$cell_line==y$cell_line,x$phenotype==y$phenotype))

# Optionally, view the result:
print(unique(most_common_virus))

common_virus_list=unique(most_common_virus)


## --------------------------------------------------------------------------------------------------------
export_list_for_plink=ind_coord_virus[,c(1,2)]
export_list_for_plink$cell_line=row.names(export_list_for_plink)
export_list_for_plink=export_list_for_plink[,c(3,1)]
#write_delim(x = export_list_for_plink,file = "plink_virus_phenotype.txt",delim = "\t")



## --------------------------------------------------------------------------------------------------------
phenotype_virus_susceptible=data.frame("IID"=toupper(row.names(susceptible_virus_clusters)),"Phenotype"=2)
phenotype_virus_resistant=data.frame("IID"=toupper(row.names(resistant_virus_clusters)),"Phenotype"=1)

phenotype_virus=rbind(phenotype_virus_susceptible,phenotype_virus_resistant)

#write_delim(x = phenotype_virus,file = "../../inputs/virus_phenotype.txt",delim = "\t")

