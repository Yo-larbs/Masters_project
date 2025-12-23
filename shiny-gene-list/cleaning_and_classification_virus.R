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
#install.packages("dplyr")
#install.packages("rstudioapi")
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(ggrepel)
library(minpack.lm)
library(drc)
library(factoextra)
library(rstudioapi)

## --------------------------------------------------------------------------------------------------------
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path))
HIV1_virus_data=read_csv('hiv1_virus_collated_n0_uncleaned(Sheet1).csv')


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
#scatter_plot_virus


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
#box_plot_virus


## --------------------------------------------------------------------------------------------------------
class(HIV1_virus_data)


## --------------------------------------------------------------------------------------------------------

HIV1_virus_data=HIV1_virus_data%>%
  group_by(screen_nb)%>%
  mutate(zscore=((((assay_output-mean(assay_output))/sd(assay_output)))))%>%
  ungroup()

virus_data_per_cell_line=HIV1_virus_data%>%nest_by(cell_line,.keep = T)





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
#logarithmic_plot_virus


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
#logistic_plot_virus



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

## --------------------------------------------------------------------------------------------------------
PCA=prcomp(x = HIV1_virus_data_PCA_1)


## ----PCA dim 1-------------------------------------------------------------------------------------------
res_var_virus=get_pca_var(PCA)
res.ind <- get_pca_ind(PCA)
dim(res.ind$coord)
ind_coord_virus <- as.data.frame(res.ind$coord)
ind_coord_virus=ind_coord_virus[-which(row.names(ind_coord_virus)==('TZM')|row.names(ind_coord_virus)=='CTR_M2_O5'|row.names(ind_coord_virus)=='CTR_M2_O5_21'),]
res.ind$cos2=res.ind$cos2[-which(row.names(res.ind$cos2)==('TZM')|row.names(res.ind$cos2)=='CTR_M2_O5'|row.names(res.ind$cos2)=='CTR_M2_O5_21'),]

