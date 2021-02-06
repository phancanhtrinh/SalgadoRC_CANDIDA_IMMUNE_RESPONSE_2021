loadPackagesCCA <- function() {
  message("\n", "Installing and loading required packages ...")

  suppressMessages({
    load_lib <- c(
    #data manipulation:
      "readxl",
      "purrr",
      "DT",
      "rlist",
      "dplyr",
    #data description:
      "MVN",
    #cca:
      "DFA.CANCOR",
      "whitening",
    #clustering:
      "textshape",
    #graphics:
      "pheatmap",
      "ggplot2",
      "ggrepel",
      "GGally",
      "sna",
      "network",
      "ggnetwork"
  )

    install_lib <- load_lib[!load_lib %in% installed.packages()]
    for (lib in install_lib)
      install.packages(lib, dependencies = TRUE)

    loading <- sapply(load_lib, require, character = TRUE)
  })

  if (!any(loading)) {
    message("\n", "Packages not working: ",  load_lib[!loading])
    stop("The above packages were not properly installed.")
  } else {
    message("\n", "Packages successfully loaded.")
  }

  return(load_lib)
}

packagesReferences <- function() {
  ref <- c(
    "whitening",  "10.1186/s12859-018-2572-9",  "",
    "ggplot2",    "",                           "331924275X"
  )

  df_ref <- as.data.frame(matrix(ref, ncol = 3, byrow = TRUE))
  colnames(df_ref) <- c("Package", "DOI", "ISBN-10")

  return(df_ref)
}

dataLoad <- function() {
  message("\n", "A FILE CHOOSER WINDOW WAS OPENED:")
  message("CHOOSE AN Excel FILE CONTAINING SPREADSHEETS WITH DATA ...")

  file <- choose.files(multi = FALSE)
  file_format <- excel_format(file)

  if (!is.na(file_format)) {
    dfs <- file %>% excel_sheets() %>% set_names() %>% map(read_excel, path = file)
    message("\n", "File successfully chosen.")
    message("\n", "File format: ", file_format, "\nNumber of spreadsheets: ", length(dfs))
  } else {
    stop("File is not an Excel file.")
  }

  output <- list(dfs = dfs, n_sheets = length(dfs))

  return(output)
}

dataLogTransform <- function(dfs, base = 2) {
  if (base <= 0) stop("Logarithm base should be positive")

  log_dfs <- lapply(dfs, function(df) {
    if (any(df < 0)) {
      stop("Data with negative entries are not allowed.")
    } else if (any(df == 0)) {
      message("\n", "Zero values identified in the dataset.")
      df <- log(df + 1, base)
      message("\n", "Data was transformed using log", base,"(value + 1)")
    } else {
      df <- log(df, base)
      message("\n", "Data was transformed using log", base,"(value)")
    }
  })
  names(log_dfs) <- names(dfs)

  return(log_dfs)
}

dataScale <- function(dfs) {
  scale_dfs <- lapply(dfs, scale)
  names(scale_dfs) <- names(dfs)
  message("\n", "Data was scaled.")

  return(scale_dfs)
}

newCANCOR <- function (data, set1, set2) {
  data <- as.data.frame(data[, c(set1, set2)])
  if (anyNA(data) == TRUE) {
    na_rows <- unique(which(is.na(data), arr.ind = TRUE)[, 1])
    data <- na.omit(data)
    warning("Incomplete cases were excluded from datasets. Missing data detected in rows: ", na_rows, ".")
  }
  else {
    message("\n", "No missing data detected in the datasets")
  }
  set1data <- as.data.frame(data[, set1])
  set2data <- as.data.frame(data[, set2])
  Ncases <- nrow(set1data)
  NVset1 <- ncol(set1data)
  NVset2 <- ncol(set2data)
  CorrelSet1 <- stats::cor(set1data)
  CorrelSet2 <- stats::cor(set2data)
  CorrelSet1n2 <- stats::cor(set2data, set1data)
  output <- DFA.CANCOR:::canonical.cor(set1data, set2data)
  cancorrels <- output$cancorrels

  CANCORoutput <- list(cancorrels = cancorrels, CoefRawSet1 = output$raw1, 
    CoefRawSet2 = output$raw2, CoefStruct11 = output$struct11, 
    CoefStruct21 = output$struct21, CoefStruct12 = output$struct12, 
    CoefStruct22 = output$struct22, CoefStandSet1 = output$stand1, 
    CoefStandSet2 = output$stand2, CorrelSet1 = CorrelSet1, CorrelSet2 = CorrelSet2, 
    CorrelSet1n2 = CorrelSet1n2, X = set1data, Y = set2data)

  return(CANCORoutput)
}

statCCA <- function(X, Y) {
  cca_data <- data.frame(X, Y)
  cca <- newCANCOR(cca_data, colnames(X), colnames(Y))

  X <- as.matrix(cca$X)
  Y <- as.matrix(cca$Y)
  x <- ncol(X)
  y <- ncol(Y)
  
  cca_cor <- whitening::cca(X, Y)$lambda

  S <- cov(cbind(X, Y))
  Sx <- S[1:x, 1:x]
  Sy <- S[(x + 1):(x + y), (x + 1):(x + y)]
  Sxy <- S[1:x, (x + 1):(x + y)]

  Rux <- cca$CoefStruct11
  Rvy <- cca$CoefStruct22

  Ru <- sum(diag(Sx))
  Rv <- sum(diag(Sy))

  if (Ru == x & Rv == y & !any(abs(signif(S, 5)) > 1)) {
    Azi2 <- lapply(as.data.frame(Rux), function(i) {
      m <- matrix(i, ncol = 1)
      m2 <- m %*% t(m)
    })
    Bzi2 <- lapply(as.data.frame(Rvy), function(i) {
      m <- matrix(i, ncol = 1)
      m2 <- m %*% t(m)
    })

    var_expl <- do.call(rbind, lapply(1:min(x, y), function(i) {
      Ru2 <- sum(diag(Reduce("+", Azi2[1:i]))) / Ru
      Rv2 <- sum(diag(Reduce("+", Bzi2[1:i]))) / Rv
      output <- c(U = Ru2, V = Rv2)
    }))
    colnames(var_expl) <- paste0("Cumulative_Variance_of_", colnames(var_expl))

    message("\n", "Consider the output 'results$var_expl' only if your data are scaled.")
  } else {
    var_expl <- NULL
  }

  output <- list(
    summary = c(n_var_X = x, n_var_Y = y),
    results = list(cca_cor = cca_cor, var_expl = var_expl),
    coefs = list(A = cca$CoefRawSet1, B = cca$CoefRawSet2),
    scores = list(U = X %*% cca$CoefRawSet1, V = Y %*% cca$CoefRawSet2),
    correlations = list(Rux = cca$CoefStruct11, Rvy = cca$CoefStruct22, 
                        Rvx = cca$CoefStruct21, Ruy = cca$CoefStruct12),
    cca_data = list(X = X, Y = Y),
    pearson_cor = list(corX = Sx, corY = Sy, corXY = Sxy)
  )

  message("\n", "CCA performed successfully.")

  return(output)
}

circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
  r <- diameter
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

newDir <- function(dir_name, at_root = TRUE, root) {
  if (at_root) root <- getwd()

  path <- paste(root, dir_name, sep = "/")

  if (!dir.exists(path)) {
    result <- dir.create(path)
    if (result) message("\n", "A folder named '", dir_name, "' was created in: ", path)
  }

  return(path)
}

saveTable <- function(p, title, ...) {
  path <- newDir(...) 

  if (dir.exists(path)) {
    write.csv(p, file = paste(path, paste0(title, ".csv"), sep = "/"), quote = FALSE, row.names = FALSE)
    DT::saveWidget(datatable(round(p, 5)), file = paste(path, paste0(title, ".html"), sep = "/"), selfcontained = FALSE)

    message("\n", "table '",  title, "' saved in: ", path, ".")
  } else {
    stop("It was not possible to create a folder to store tables.")
  }
}

savePlot <- function(p = NULL, title, width, height, extension = "pdf", ...) {
  path <- newDir(...) 

  if (dir.exists(path) & !is.null(p)) {
    suppressWarnings({
      if(extension == "pdf") {
        pdf_name <- paste0(title, ".pdf")
        pdf(paste(path, pdf_name, sep = "/"), width = width, height = height)
        print(p)
        dev.off()
      } else {
        png_name <- paste0(title, ".png")
        png(paste(path, png_name, sep = "/"), width = width, height = height, units = "in", res = 300)
        print(p)
        dev.off()
      }
    }) 

    message("\n", "plot '", title , "' saved in: ", path, ".")
  } else if (is.null(p)) {
    return(paste(path, paste(title, extension, sep = "."), sep = "/"))
  } else {
    stop("It was not possible to create a folder to store figures.")
  }
}

plotCCA <- function(cca, canonical_variate = 1, k_min = 0.6, x_title = "X", y_title = "Y", short_x_title = "X",
  short_y_title = "Y", x_col = "red", y_col = "blue") {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))
  k_min <- abs(k_min)

  if ((canonical_variate < 0 | canonical_variate > (max_cvs - 1))) stop("canonical_variate should be at least equal to 1 and most equal to ", max_cvs - 1, ".")
  if (k_min > 1) stop("k_min should be at between 0 and 1, inclusive.")

  cca_x <- cca$correlations$Rux[, canonical_variate:(canonical_variate + 1)]
  df_cca_x <- data.frame(cca_x, Vars = colnames(cca$cca_data$X), group = "X",
   col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_x[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 1))
   }), check.names = FALSE)
  colnames(df_cca_x)[1:2] <- paste0("CV", 1:2)

  cca_y <- cca$correlations$Rvy[, canonical_variate:(canonical_variate + 1)]
  df_cca_y <- data.frame(cca_y, Vars = colnames(cca$cca_data$Y), group = "Y",
   col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_y[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 2))
   }))
  colnames(df_cca_y)[1:2] <- paste0("CV", 1:2)

  df_cca <- rbind(df_cca_x, df_cca_y)

  facet_label <- c(X = x_title, Y = y_title)

  df_cca$label_x <- df_cca$label_y <- NA
  df_cca$label_x[grep("X", df_cca$group)[1]] <- paste0(short_x_title, "-CV", canonical_variate)
  df_cca$label_y[grep("X", df_cca$group)[1]] <- paste0(short_x_title, "-CV", canonical_variate + 1)
  df_cca$label_x[grep("Y", df_cca$group)[1]] <- paste0(short_y_title, "-CV", canonical_variate)
  df_cca$label_y[grep("Y", df_cca$group)[1]] <- paste0(short_y_title, "-CV", canonical_variate + 1)

  p <- ggplot(data = df_cca, aes(x = CV1, y = CV2, colour = col, group = group)) +
   geom_point(alpha = 1, size = 2, show.legend = FALSE) +
   geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = 1, colour = "black", fill = "transparent",
    linetype = "dotted") +
   geom_rect(xmin = -k_min, xmax = k_min, ymin = -k_min, ymax = k_min, colour = "black",
    fill = "transparent", linetype = "dotted") +
   geom_hline(yintercept = 0, alpha = 0.6, color = "gold", size = 1.5) +
   geom_vline(xintercept = 0, alpha = 0.6, color = "gold", size = 1.5) +
   geom_text(aes(label = label_x), x = 0.2, y = 0.05, color = "gold", size = 4.2, fontface = "bold") +
   geom_text(aes(label = label_y), x = -0.05, y = -0.2, color = "gold", angle = 90, size = 4.2,
    fontface = "bold") +
   facet_wrap(~group, labeller = as_labeller(facet_label)) +
   coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
   scale_colour_manual(values = c("0" = "grey", "1" = x_col, "2" = y_col)) +
   geom_text_repel(aes(label = ifelse(col != "0", Vars, NA)), size = 3.3,
   segment.alpha = 0.4, segment.size = 0.3, min.segment.length = 0, fontface = "bold",
    show.legend = FALSE) +
   theme_classic() +
   theme(
    legend.text = element_text(size = 11),
    axis.line = element_blank(),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 13)) +
   labs(y = paste("Correlation with canonical variate", canonical_variate + 1),
    x = paste("Correlation with canonical variate", canonical_variate))

  output <- list(p = p, df = df_cca)

  return(output)
}

plotCorCCA <- function(cca, n = 2, k = 0.86, seed = 100, node_size = 5, label_size = 2.3, x_col = "red", y_col = "blue") {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))

  if (n < 0 | n > max_cvs) stop("n should be at least equal to 1 and at most equal to ", max_cvs, ".")
  if (k > 1) stop("k should be at between 0 and 1, inclusive.")

  X <- cca$cca_data$X
  Y <- cca$cca_data$Y
  signed_cor <- cca$results$cca_cor

  cor_Vars_cca <- cor(cbind(X, Y))

  n <- min(n, length(signed_cor))
  n_cor <- seq_len(n)

  cor_cv <- signed_cor[n_cor]

  u_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_ux_cv <- cca$correlations$Rux[, i]
    cor_uy_cv <- cca$correlations$Rvx[, i]
    u_cv <- c(cor_ux_cv, cor_uy_cv)
  }))

  v_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_vx_cv <- cca$correlations$Ruy[, i]
    cor_vy_cv <- cca$correlations$Rvy[, i]
    v_cv <- c(cor_vx_cv, cor_vy_cv)
  }))

  diag_cor <- diag(cor_cv)
  diag_1 <- diag(2 * n)
  diag_1[1:n, (n + 1):ncol(diag_1)] <- diag_cor
  diag_1[(n + 1):nrow(diag_1), 1:n] <- diag_cor

  cv_cols <- cbind(u_cv, v_cv)
  cv_rows <- cbind(t(cv_cols), diag_1)
  bind_cols <- cbind(cor_Vars_cca, cv_cols)
  bind_rows <- rbind(bind_cols, cv_rows)

  colnames(bind_rows) <- rownames(bind_rows) <- c(
   colnames(cor_Vars_cca), paste0("X-CV", n_cor), paste0("Y-CV", n_cor)
  )

  diag(bind_rows) <- 0

  message("\n", "correlation matrix is symmetric? ", isSymmetric(bind_rows))

  cor_Vars_cv_k <- bind_rows
  cor_Vars_cv_k[abs(cor_Vars_cv_k) < k] <- 0
  cor_Vars_cv_k[abs(cor_Vars_cv_k) >= k] <- 1

  type_Vars <- ifelse(
   colnames(bind_rows) %in% colnames(X), "X",
   ifelse(colnames(bind_rows) %in% colnames(Y), "Y", "CV")
  )
  name_Vars <- colnames(bind_rows)
  name_Vars[colSums(cor_Vars_cv_k) == 0] <- ""

  net <- network(cor_Vars_cv_k, directed = FALSE)
  net %v% "type" <- type_Vars
  net %v% "name" <- name_Vars

  set.seed(seed)
  df_network <- ggnetwork(net, layout = "kamadakawai")

  p <- ggplot(data = df_network, aes(x, y, xend = xend, yend = yend)) +
   geom_edges(alpha = 0.8, color = "grey") +
   geom_nodes(aes(color = type), alpha = 0.5, size = node_size, show.legend = FALSE) +
   geom_nodetext(aes(label = name), fontface = "bold", size = label_size, alpha = 0.7) +
   scale_colour_manual("", values = c("X" = x_col, "Y" = y_col, "CV" = "gold")) +
   theme_blank()

  output <- list(p = p, df_network = df_network)

  return(output)
}

checkParameters <- function(list_param) {
  parameters <- list(
    empty = FALSE,
    heatmap_clust_method = "ward.D2",
    heatmap_color_ramp = rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")),
    heatmap_border_color = grey(0.4),
    heatmap_fontsize = 5,
    plotCCA_canonical_variate = 1,
    plotCCA_corr_min = 0.8,
    plotCCA_x_title = "Dataset X",
    plotCCA_y_title = "Dataset Y",
    plotCCA_short_x_title = "U",
    plotCCA_short_y_title = "V",
    plotCorCCA_n_canonical_variates = 2,
    plotCorCCA_corr_min = 0.8,
    plotCorCCA_seed = 1,
    plotCorCCA_node_size = 5,
    plotCorCCA_label_size = 3,
    plots_x_color = "red",
    plots_y_color = "blue",
    plots_file_extension = "pdf",
    plots_width = 13,
    plots_height = 7,
    data_log_trans = FALSE,
    data_log_base = 2,
    data_scale_trans = TRUE,
    data_sheet1 = 1,
    data_sheet2 = 2,
    data_X = NA,
    data_Y = NA
  )

  if (!is.null(list_param)) {
    user_param <- names(list_param)
    cond <- any(is.na(match(user_param, names(parameters))))
    if (cond) {
      stop("there are some undefined parameters on the list provided through the argument in list_param.")
    } else {
      for (i in user_param) parameters[[i]] <- list_param[[i]]
    }
  }

  return(parameters)
}

runCCA <- function(list_param = list(empty = TRUE), publish = FALSE) {
  parameters <- checkParameters(list_param)
  for (i in names(parameters)) assign(i, parameters[[i]])

  loaded_packages <- loadPackagesCCA()

  new_path <- newDir("output_CCA")

  if(!publish) {
    current_output <- gsub(":| ", "-", format(Sys.time(), "%Y-%m-%d_%X"))
    new_path <- newDir(current_output, at_root = FALSE, root = new_path)

    list.save(parameters, file = paste(new_path, "parameters.json", sep = "/"))
    cond0 <- "parameters.json" %in% list.files(new_path, pattern = "^(parameters.json)$")
    if(!cond0) stop("Could not save the file 'parameters.json' in ", new_path, ".")
  }

  if (publish) {
    new_path <- newDir("R_CCA", at_root = FALSE, root = new_path)
  }

  if (!is.data.frame(data_X) | !is.data.frame(data_Y)){
    message("\n", "One of the datasets, X or Y, are not data.frame or were set as NA.")

    dfs <- dataLoad()
    n_sheets <- dfs$n_sheets

    if (n_sheets < 2) {
      stop("Excel file should have at least 2 spreadsheets.")
    } else if ((data_sheet1 > 0 & data_sheet1 <= n_sheets) & (data_sheet2 > 0 & data_sheet2 <= n_sheets)) {
      message("Spreadsheets ", data_sheet1, " and ", data_sheet2, " were selected for CCA.")
    } else {
      data_sheet1 = 1
      data_sheet2 = 2
      warning("reset variables to: data_sheet1 = 1 and data_sheet2 = 2.")
      warning("Only the first two spreadsheets were used for CCA.")
    }

    dfs$dfs <- dfs$dfs[c(data_sheet1, data_sheet2)]
  } else {
    dfs <- list(dfs = list(data_X, data_Y))
    n_sheets <- NULL
  }

  dimXY <- sapply(dfs$dfs, dim)
  if (diff(dimXY[1, ])) stop("Datasets does not have the same number of rows.")
  if (any(dimXY[2, ] == 0)) stop("Each dataset must have at least one variable.")

  if (data_log_trans) dfs$dfs <- dataLogTransform(dfs$dfs, base = data_log_base)
  if (data_scale_trans) dfs$dfs <- dataScale(dfs$dfs)

  data_X <- dfs$dfs[[1]]
  data_Y <- dfs$dfs[[2]]

  cca <- statCCA(data_X, data_Y)

  cca_results <- round(data.frame(cca$results$var_expl, Canonical_Correlation = cca$results$cca_cor), 5)
  saveTable(cca_results, "cca_results", dir_name = "tables", at_root = FALSE, root = new_path)

  if (publish) {
    saveTable(cca$cca_data$X, "dataset_X", dir_name = "data", at_root = FALSE, root = new_path)
    saveTable(cca$cca_data$Y, "dataset_Y", dir_name = "data", at_root = FALSE, root = new_path)

    path_code <- newDir("code", at_root = FALSE, root = new_path)

    cond1 <- "CCApackage.R" %in% list.files(path_code, pattern = "^(CCApackage.R)$")
    if (!cond1) {
      cond2 <- length(list.files(getwd(), pattern = "^(CCApackage.R)$"))
      if (!cond2) stop("file named 'CCApackage.R' not found.")

      file_to <- paste(path_code, "CCApackage.R", sep = "/")
      file_from <- paste(getwd(), "CCApackage.R", sep = "/")

      file.create(file_to)
      filecopy <- file.copy(from = file_from, to = file_to, overwrite = TRUE)
      if (!filecopy) stop("Could not copy file 'CCApackage.R' from workspace to ", path_code, ".")
    }

    parameters[["data_log_trans"]] <- FALSE
    parameters[["data_scale_trans"]] <- FALSE
    parameters[["data_X"]] <- as.data.frame(cca$cca_data$X)
    parameters[["data_Y"]] <- as.data.frame(cca$cca_data$Y)
    list.save(parameters, file = paste(path_code, "list_param.rdata", sep = "/"))
    cond3 <- length(list.files(path_code, pattern = "^(list_param.rdata)$"))
    if(!cond3) stop("Could not save file 'list_param.rdata' in ", path_code, ".")

    write_package <- "if (!require('rlist')) install.packages('rlist')\nlibrary(rlist)"
    write_path <- "path <- paste(getwd(), 'code', sep = '/')"
    write_source <- "source(paste(path, 'CCApackage.R', sep = '/'))"
    write_param <- "list_param <- list.load(paste(path, 'list_param.rdata', sep = '/'))"
    write_runCCA <- "runCCA(list_param, publish = FALSE)"
    write.table(c(write_package, write_path, write_source, write_param, write_runCCA), file = paste(path_code, "main.R", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    cond4 <- length(list.files(path_code, pattern = "^(main.R)$"))
    if(!cond4) stop("Could not save file 'main.R' in ", path_code, ".")

    path_ref <- newDir("references-and-packages", at_root = FALSE, root = new_path)
    write.csv(packagesReferences(), file = paste(path_ref, "main_references.csv", sep = "/"), quote = FALSE, row.names = FALSE)
    write.table(data.frame(loaded_packages), file = paste(path_ref, "used_packages.csv", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    readme <-  "
      1 - Start a session in R version > 4.0.1\n
      2 - Difene the folder 'R_CCA' as your working directory\n
      3 - Open the file 'main.R' located in the folder 'code'\n
      4 - Excute all commands in 'main.R'\n
      The results will be saved in a folder named 'output_CCA'."
    write.table(readme, file = paste(new_path, "READ-ME.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  cor_clust <- cluster_matrix(cca$pearson_cor$corXY, dim = "both", method = heatmap_clust_method)
  p1 <- pheatmap(cor_clust, cluster_rows = FALSE, cluster_cols = FALSE, 
    color = heatmap_color_ramp, border_color = heatmap_border_color, fontsize = heatmap_fontsize,
    filename = savePlot(title = "heatmap_of_X_against_Y", extension = plots_file_extension, dir_name = "figures", at_root = FALSE, root = new_path)
  )
  
  p2 <- plotCCA(cca, plotCCA_canonical_variate, plotCCA_corr_min, plotCCA_x_title, plotCCA_y_title, plotCCA_short_x_title, plotCCA_short_y_title, plots_x_color, plots_y_color)
  savePlot(p2$p, "correlations-of-variables-against-canonical-variables", plots_width, plots_height, plots_file_extension, dir_name = "figures", at_root = FALSE, root = new_path)

  p3 <- plotCorCCA(cca, plotCorCCA_n_canonical_variates, plotCorCCA_corr_min, plotCorCCA_seed, plotCorCCA_node_size, plotCorCCA_label_size, plots_x_color, plots_y_color)
  savePlot(p3$p, "graph-of-pearson-correlations-including-canonical-variables", plots_width, plots_height, plots_file_extension, dir_name = "figures", at_root = FALSE, root = new_path)

  if (!publish) {
    mvn_X <- mvn(cca$cca_data$X)
    mvn_Y <- mvn(cca$cca_data$Y)

    desc_X <- arrange(round(mvn_X$Descriptives, 5), Skew, Kurtosis)
    desc_Y <- arrange(round(mvn_Y$Descriptives, 5), Skew, Kurtosis)
    n_var <- floor(nrow(cca$cca_data$X) / 20)

    mvn_X$Descriptives <- desc_X
    mvn_Y$Descriptives <- desc_Y
    
    saveTable(desc_X, "description_X", dir_name = "descriptive-analysis", at_root = FALSE, root = new_path)
    saveTable(desc_Y, "description_Y", dir_name = "descriptive-analysis", at_root = FALSE, root = new_path)

    cat("Dataset X -------------------------------------------------------------------------", "\n\n")
    print(mvn_X)
    cat("Dataset Y -------------------------------------------------------------------------", "\n\n")
    print(mvn_Y)
    cat("Suggested Number of Variables -----------------------------------------------------", "\n\n")
    cat("It is suggested to use at most ", n_var, " variables in this CCA.", "\n")
    cat("User has chosen to keep ", nrow(desc_X) + nrow(desc_Y), " variables.", "\n\n")
    cat("-----------------------------------------------------------------------------------", "\n")
  }
}

