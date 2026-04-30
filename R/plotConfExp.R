#' Visualize confounder level distribution by (binary) exposure level
#'
#' @inheritParams devMSMHelpers_common_docs
#' 
#' @param confounders list of names of time-varying and time-invariant confounders that 
#' contain at least one integer/factor 
#' @param data single complete wide dataset
#' @param exposures list of exposure variable names
#' @param save.path file name path if you want to save the figure locally
#' @return none
#' @export
#' @examples
#' library(devMSMs)
#' test <- data.frame(ID = 1:50,
#'                    A.1 = rnorm(n = 50),
#'                    A.2 = rnorm(n = 50),
#'                    A.3 = rnorm(n = 50),
#'                    B.1 = rnorm(n = 50),
#'                    B.2 = rnorm(n = 50),
#'                    B.3 = rnorm(n = 50),
#'                    C = rnorm(n = 50),
#'                    D.3 = rnorm(n = 50))
#' obj <- initMSM(
#'   test,
#'   exposure = c("A.1", "A.2", "A.3"),
#'   ti_conf = c("C"),
#'   tv_conf = c("B.1", "B.2", "B.3", "D.3")
#' )
#'
#'plotConfExp(confounders = c(ti_conf, tv_conf), data = test, exposures = c("A.1", "A.2", "A.3"))


plotConfExp <- function(confounders, data, exposures, save.path = NULL){
  
  dreamerr::check_arg(data, "data.frame")
  dreamerr::check_arg(confounders, "vector character")
  dreamerr::check_arg(exposures, "vector character")
  dreamerr::check_arg(save.path, "vector character | NULL")
  
  
  # find factor + integer confounders
  factor_integer_confounders <- confounders[
    sapply(data[confounders], function(x) is.factor(x) || is.integer(x) || is.factor(x))
  ]
  
  if(is.null(length(factor_integer_confounders))){
    stop("There are no integer confounders to plot.")
  }
  
  plot_data <- data %>%
    tidyr::pivot_longer(
      cols = all_of(exposures),
      names_to = "exposure",
      values_to = "exposure_level",
      values_transform = list(exposure_level = as.character)
    ) %>%
    dplyr::mutate(exposure_level = factor(exposure_level)) %>%
    tidyr::pivot_longer(
      cols = all_of(factor_integer_confounders),
      names_to = "confounder",
      values_to = "level",
      values_transform = list(level = as.character)
    ) %>%
    dplyr::count(exposure, confounder, exposure_level, level)
  
  plot_data <- plot_data %>%
    dplyr::group_by(confounder, exposure) %>%   # match your facet structure
    dplyr::mutate(
      max_n = max(n, na.rm = TRUE),
      label_x = n + 0.05 * max_n        # 5% padding
    ) %>%
    ungroup()
  
  # compute smallest cell count per confounder (across all exposures)
  min_cell_df <- plot_data %>%
    dplyr::group_by(confounder) %>%
    dplyr::summarise(min_cell = min(n, na.rm = TRUE), .groups = "drop")
  
  plot_data <- plot_data %>%
    dplyr::left_join(min_cell_df, by = "confounder") %>%
    dplyr::mutate(
      confounder = fct_reorder(confounder, min_cell),
      exposure = factor(exposure, levels = exposures)
    )
  
  dodge <- position_dodge(width = 0.9)
  
  p <- ggplot(plot_data, aes(y = level, x = n, fill = exposure_level, group = exposure_level)) +
    geom_col(position = dodge) +
    geom_text(
      aes(
        x = label_x,
        label = ifelse(n < 20, n, "")
      ),
      position = dodge,
      hjust = 0,
      size = 3,
      color = "red"
    ) +
    scale_fill_manual(
      values = c("0" = "#6B8E23",
                 "1" = "#C97C5D")
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    facet_grid(confounder ~ exposure, scales = "free_y") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      strip.text.y = element_text(angle = 0)
    ) +
    labs(
      x = "Count",
      y = "Confounder Level",
      fill = "Exposure Level",
      title = "Counts of each Factor/Integer Confounder Levels at Each Exposure Level"
    )
  if(!is.null(save.path)){
    ggsave(file.path(save.path), width = 12, height = 10)
  }
  p
}
