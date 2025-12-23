## ----setup, include=FALSE--------------------------------------------------------------------------------
knitr::opts_chunk$set(
  echo = TRUE,
  dev = ifelse(knitr::is_latex_output(), "pdf", "svg"),
  fig.ext = ifelse(knitr::is_latex_output(), "pdf", "svg"),
  fig.path = "../../output/figures/cleaning_and_classification/vector/"
  # write files with .svg extension
  )



## --------------------------------------------------------------------------------------------------------
#install.packages("factoextra")
#install.packages("tidyr")
#install.packages("readr")
#install.packages("ggplot2")
#install.packages("DescTools")
#install.packages('ggrepel')
#install.packages("minpack.lm")
#install.packages("drc")
#install.packages("purrr")
library(tidyr)
library(readr)
library(dplyr)
library(ggplot2)
library(DescTools)
library(ggrepel)
library(minpack.lm)
library(drc)
library(factoextra)
library(purrr)
library(rstudioapi)


## --------------------------------------------------------------------------------------------------------
current_path = rstudioapi::getActiveDocumentContext()$path
setwd(dirname(current_path ))
HIV1_vector_data=read_csv('HIV1_vectors_collated_n0_uncleaned(Sheet1).csv')


## --------------------------------------------------------------------------------------------------------
HIV1_vector_data=HIV1_vector_data[,c(2,5,6,7:13,15)]

HIV1_vector_data=na.omit(HIV1_vector_data)
HIV1_vector_data=HIV1_vector_data%>%mutate(cell_line=if_else(cell_line=="Blix_1","Bilx_1",HIV1_vector_data$cell_line))


## --------------------------------------------------------------------------------------------------------
library(dplyr)
HIV1_vector_data=HIV1_vector_data%>%group_by(cell_line,screen_nb,replicate)%>% 
  mutate(max_gfp=max(assay_output),max_gfp_titre=which.max(assay_output)) %>%
mutate(assay_output=if_else((infection_volume_ul>infection_volume_ul[max_gfp_titre] & assay_output  <= 0.8 * max_gfp),max_gfp,assay_output))%>% dplyr::select(-max_gfp,-max_gfp_titre) %>% ungroup()


## --------------------------------------------------------------------------------------------------------
vector_data_per_cell_line=HIV1_vector_data%>%nest_by(cell_line,.keep = T)


## --------------------------------------------------------------------------------------------------------
apply_plot=function(i){ggplot(data = i,aes(x=infection_volume_ul,y=assay_output,group = interaction(screen,replicate) ))+
  geom_point()+
  ggtitle(i$cell_line[1])+
  theme_classic()
}


## ----scatter_plot----------------------------------------------------------------------------------------
scatter_plot_vector=purrr::map(vector_data_per_cell_line$data,apply_plot)
names(scatter_plot_vector)=vector_data_per_cell_line$cell_line
#scatter_plot_vector


## --------------------------------------------------------------------------------------------------------
apply_plot_boxplot=function(i){ggplot(data = i,aes(y=assay_output,x=titre,))+
  geom_boxplot()+
  ggtitle(i$cell_line[1])+
  #stat_summary(geom="text", fun.y=quantile,
               #aes(label=sprintf("%1.1f", ..y..), color=factor(titre)),
               #position=position_nudge(x=0.6), size=3.5) +
  theme_classic()
}


## ----box_plot--------------------------------------------------------------------------------------------
box_plot_vector=purrr::map(vector_data_per_cell_line$data,apply_plot_boxplot)
names(box_plot_vector)=vector_data_per_cell_line$cell_line
#box_plot_vector


## ----unnormalised_boxplot--------------------------------------------------------------------------------
unormalised_boxplot=ggplot(data = HIV1_vector_data,aes(y=assay_output, x=as.factor(screen_nb),group = screen_nb,))+
  geom_boxplot()+
  ggtitle('change between screens')+
  theme_classic()
#unormalised_boxplot


## --------------------------------------------------------------------------------------------------------

HIV1_vector_data=HIV1_vector_data%>%
  group_by(screen_nb)%>%
  mutate(zscore=((((assay_output-mean(assay_output))/sd(assay_output)))))%>%
  ungroup()

vector_data_per_cell_line=HIV1_vector_data%>%nest_by(cell_line,.keep = T)



## ----normalised_boxplot----------------------------------------------------------------------------------
normalised_boxplot=ggplot(data = HIV1_vector_data,aes(y=zscore, x=as.factor(screen_nb),group = screen_nb,))+
  geom_boxplot()+
  ggtitle('change between screens')+
  theme_classic()
#normalised_boxplot


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
all_outliers=purrr::map(vector_data_per_cell_line$data,outlier_function)

HIV1_vector_data_outliers=bind_rows(all_outliers)%>%
  group_by(cell_line,screen_nb,replicate)%>%
  filter(((sum(is.na(outlier)))<3))

HIV1_vector_data=bind_rows(all_outliers)%>%
  group_by(cell_line,screen_nb,replicate)%>%
  filter(((sum(is.na(outlier)))>2))




## --------------------------------------------------------------------------------------------------------
HIV1_vector_data=HIV1_vector_data%>%group_by(cell_line,screen,screen_nb,replicate)%>%
  mutate(area_under_curve=AUC(infection_volume_ul,(zscore)))%>%
  ungroup()


## --------------------------------------------------------------------------------------------------------
vector_data_per_cell_line=HIV1_vector_data%>%nest_by(cell_line,.keep = T)


## ----message=FALSE,include=FALSE-------------------------------------------------------------------------
logarithmic_func <- function(x, a, b) {
  a + b * log(x,base = 2)
}
fit_log <- function(data) {

  # Fit the logarithmic model for each group and add coefficients to the data frame
  fitted_data <- data %>%
    group_by(screen_nb,replicate) %>%
    group_modify(~ {
      tryCatch({
        # Fit the logarithmic model
        model <- nlsLM(
          zscore ~ logarithmic_func(infection_volume_ul, a, b),
          data = .x,
          start = list(a = 0, b = -1),# c = -0.0001),
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
vector_data_per_cell_line=purrr::map(vector_data_per_cell_line$data,fit_log)%>%bind_rows()%>%nest_by(cell_line,.keep = T)


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


## ----logarithmic_plot------------------------------------------------------------------------------------
logarithmic_plot_vector=purrr::map(vector_data_per_cell_line$data,apply_plot_log)
names(logarithmic_plot_vector)=vector_data_per_cell_line$cell_line
#logarithmic_plot_vector


## ----include=FALSE---------------------------------------------------------------------------------------
# Correct 4PL formula (parentheses fix)
logistic_func <- function(x, b, d, e,c){
 c+ ((d-c)/(1 + exp(b*(x -e))))# Fixed denominator placement
}

# Robust fitting function
logistic_fit <- function(data) {
    fitted_data <- data %>%
        group_by(screen_nb,replicate) %>%  # Verify these columns exist in your data
        group_modify(~ {
            tryCatch({ # Fit with drm, explicitly stating NO log transformation
                model <- suppressWarnings(drm(zscore ~ infection_volume_ul, fct = L.4(),lowerl = c(-Inf, min(.x$zscore), -Inf, -Inf), data = .x,control = drmc(method = "CG")))
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
vector_data_per_cell_line <- purrr::map(vector_data_per_cell_line$data,logistic_fit)%>%bind_rows()%>%nest_by(cell_line,.keep = T)




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
      aes(color = factor(screen_nb)),
      color = "blue",# Color by screen_nb
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
      x = "Infection Volume (microL)",
      y = "Mean Z-Score"
    ) +
    theme_classic(base_size = 12)
  
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


## ----logistic_plot---------------------------------------------------------------------------------------
logistic_plot_vector=purrr::map(vector_data_per_cell_line$data,apply_plot_logistic)
names(logistic_plot_vector)=vector_data_per_cell_line$cell_line
#logistic_plot_vector



## --------------------------------------------------------------------------------------------------------
HIV1_vector_data=vector_data_per_cell_line%>%unnest()

HIV1_vector_data_mean=HIV1_vector_data%>%
  group_by(batch,cell_line,screen,screen_nb,titre,infection_volume_ul)%>%
  reframe(assay_output=mean(assay_output),
            zscore=mean(zscore),
            a_log=mean(a),
            b_log=mean(b),
            #c_log=mean(c),
            logis_b=mean(logis_b),
            logis_d=mean(logis_d),
            logis_c=mean(logis_c),
            logis_e=mean(logis_e),
            area_under_curve=mean(area_under_curve))

HIV1_vector_data_mean <- HIV1_vector_data_mean %>%
  mutate(across(c(zscore, a_log, b_log,#c_log, 
                  logis_b, logis_d,logis_c,logis_e, area_under_curve), as.numeric))



## --------------------------------------------------------------------------------------------------------
HIV1_vector_data_PCA=unique(HIV1_vector_data_mean%>%
                              group_by(cell_line,screen_nb)%>%
                              reframe(assay_output=max(zscore),
                                        a_log=a_log,
                                        b_log=(b_log),
                                        #c_log=c_log,
                                        logis_b=logis_b,
                                        logis_d=logis_d,
                                        logis_c=logis_c,
                                        logis_e=logis_e,
                                        area_under_curve=area_under_curve))


## --------------------------------------------------------------------------------------------------------

HIV1_vector_data_PCA_sum <- HIV1_vector_data_PCA %>%
  group_by(cell_line) %>%
  reframe(
    assay_output      = mean(assay_output, na.rm = TRUE),
    a_log             = mean(a_log, na.rm = TRUE),
    b_log             = mean(b_log, na.rm = TRUE),
    #c_log            = mean(c_log, na.rm = TRUE),
    logis_b           = -mean(logis_b, na.rm = TRUE),
    logis_d           = mean(logis_d, na.rm = TRUE),
    #logis_c           = mean(logis_c, na.rm = TRUE),
    #logis_e           = mean(logis_e, na.rm = TRUE),
    area_under_curve  = abs(mean(area_under_curve, na.rm = TRUE))
  )


## --------------------------------------------------------------------------------------------------------
HIV1_vector_data_PCA_1=HIV1_vector_data_PCA_sum%>%mutate_if(is.numeric,scale)


HIV1_vector_data_PCA_1=subset(HIV1_vector_data_PCA_1, select=-cell_line)

HIV1_vector_data_PCA_1=as.matrix(HIV1_vector_data_PCA_1)

row.names(HIV1_vector_data_PCA_1)=HIV1_vector_data_PCA_sum$cell_line





## --------------------------------------------------------------------------------------------------------
HIV1_vector_data_PCA_1=as.data.frame(HIV1_vector_data_PCA_1)
#cor(x=HIV1_vector_data_PCA_1$logis_c,y=HIV1_vector_data_PCA_1$assay_output)



## --------------------------------------------------------------------------------------------------------
PCA=prcomp(x = HIV1_vector_data_PCA_1)




## ----PCA dim 1-------------------------------------------------------------------------------------------
res_var=get_pca_var(PCA)
res.ind <- get_pca_ind(PCA)
dim(res.ind$coord)
ind_coord <- as.data.frame(res.ind$coord)
ind_coord=ind_coord[-which(row.names(ind_coord)==('293T')|row.names(ind_coord)=='CTR_M2_O5'|row.names(ind_coord)=='CTR_M2_O5_21'),]
res.ind$cos2=res.ind$cos2[-which(row.names(res.ind$cos2)==('293T')|row.names(res.ind$cos2)=='CTR_M2_O5'|row.names(res.ind$cos2)=='CTR_M2_O5_21'),]

