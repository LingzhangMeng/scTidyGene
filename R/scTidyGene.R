
#' Generate a tidy list of conserved expressed genes
#' @param seu_obj A Seurat object containing single-cell RNA-seq data.
#' @param ident.1 The identity of the cluster to analyze.
#' @param grouping.var The variable to group by when finding conserved markers.
#'
#' @return A data frame of conserved markers across clusters.
#' @export
#'
#' @examples
#' # Example usage will be added in a future release.

#' @importFrom Seurat FindConservedMarkers
#' @importFrom dplyr mutate, select, bind_rows


scCEGs <- function(seu_obj, ident.1, grouping.var) {

# Get the unique clusters
clusters <- unique(levels(seu_obj@active.ident))

# Initialize a list to store results
conserved_markers_list <- list()

# Loop through each cluster and find conserved markers
for (cluster in clusters) {

  # Find conserved markers for the current cluster
  markers <- Seurat::FindConservedMarkers(seu_obj, ident.1 = cluster, grouping.var = grouping.var, verbose = T)

  # Add Cell_Cluster and gene in the first 2 columns
  markers <- markers %>% dplyr::mutate(Cell_Cluster = cluster, gene = rownames(markers)) %>%
    dplyr::select(Cell_Cluster, gene, everything())

  # remove rownames
  rownames(markers) <- NULL

  # Add the results to the list with the cluster name
  conserved_markers_list[[as.character(cluster)]] <- markers
}

# Standardize column names and fill missing columns with NA
max_col_names <- unique(unlist(lapply(conserved_markers_list, colnames)))

# Create a function to ensure all data frames have the same columns
standardize_columns <- function(df, col_names) {
  missing_cols <- setdiff(col_names, colnames(df))
  for (col in missing_cols) {
    df[[col]] <- NA  # Add missing columns with NA
  }
  return(df[col_names])  # Return the data frame with standardized columns
}

# Apply the function to standardize columns for each data frame in the list
standardized_list <- lapply(conserved_markers_list, standardize_columns, max_col_names)

# Combine all results into a single data frame
combined_results <- dplyr::bind_rows(standardized_list)

return(combined_results)
}




#' Generate a tidy list of differentially expressed genes
#' @param seu_obj A Seurat object containing single-cell RNA-seq data.
#' @param ident.1 The identity of the first condition to compare.
#' @param ident.2 The identity of the second condition to compare.
#' @param group.by The variable to group by when performing differential expression analysis.
#' @param only.pos Logical, if TRUE, only positive markers are returned.
#' @param min.pct Minimum percentage of cells expressing the gene.
#' @param logfc.threshold Minimum log-fold change threshold.
#'
#' @return A data frame of differentially expressed genes across clusters.
#' @export
#' @examples
#' # Example usage will be added in a future release.
#' # result <- scDEGs(seu_obj, "ConditionA", "ConditionB", group.by = "group_var")

#' @importFrom Seurat FindMarkers
#' @importFrom dplyr bind_rows
#' @importFrom progress progress_bar
#' @importFrom dplyr relocate

scDEGs <- function(seu_obj, ident.1, ident.2, group.by, only.pos = TRUE,
                   min.pct = 0.25, logfc.threshold = 0.25) {

  # Validate input
  if (!inherits(seu_obj, "Seurat")) {
    stop("seu_obj must be a Seurat object.")
  }

  # Initialize an empty list to store markers for each cluster
  cluster_markers <- list()

  # Get the unique clusters
  clusters <- unique(levels(seu_obj@active.ident))

  # Initialize the progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent in :elapsed",
    total = length(clusters),
    clear = FALSE,
    width = 60
  )

  # Loop through each cluster
  for (cluster in clusters) {
    message(paste("Calculating Differentially Expressed Genes (DEG) in cluster -", cluster, "\n"))
    pb$tick()

    # Subset the Seurat object to include only cells from the current cluster
    cluster_cells <- subset(seu_obj, idents = cluster)

    # Check cell counts
    cell_counts <- as.data.frame(table(cluster_cells@meta.data[[group.by]]))
    min.N <- min(cell_counts$Freq)

    # Check if both conditions have at least 3 cells
    if (nrow(cell_counts) >= 2 && min.N >= 3) {
      # Perform differential expression analysis
      markers <- Seurat::FindMarkers(cluster_cells, ident.1 = ident.1, ident.2 = ident.2,
                                     group.by = group.by,
                                     only.pos = only.pos,
                                     min.pct = min.pct,
                                     logfc.threshold = logfc.threshold)

      # Add columns for cluster identifier and gene names
      markers$Cell_cluster <- as.character(cluster)

      markers$gene <- rownames(markers)

      rownames(markers) <- NULL

      # Store the results in the list
      cluster_markers[[as.character(cluster)]] <- markers

      message(paste("Completed processing cluster -", cluster, "\n"))
    } else {
      message(paste("Skipping cluster -", cluster, "due to insufficient cell numbers in one or both conditions.", "\n"))
    }
  }

  # Combine all cluster markers into a single data frame
  combined_markers <- dplyr::bind_rows(cluster_markers)

  rownames(combined_markers) <- NULL

  # Re-arrange columns: move Cell_clust to 1st column, gene to 2ed column
  combined_markers <- combined_markers %>% dplyr::relocate(Cell_cluster, .before = 1) %>% dplyr::relocate(gene, .after = Cell_cluster)

  return(combined_markers)
}
