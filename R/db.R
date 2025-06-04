#' Query taxonomy fields for a list of MSPs from the GMT database
#'
#' Retrieves taxonomic information for one or more MSPs (Metagenomic Species Pan-genomes)
#' from the GMT database API, with optional selection of specific fields.
#'
#' @param msp_names A character vector of MSP identifiers (e.g., \code{"msp_0005"}).
#' @param gtdb_version A numeric indicating the GTDB version to use (default: \code{220.0}). 
#' Must be one of \code{214.0} or \code{220.0}.
#' @param fields A character vector of field names to retrieve (default: empty vector, returns all fields).
#'
#' @return A named list where each element contains the taxonomy information for the corresponding MSP.
#' @export
#'
#' @examples
#' get_taxonomy(msp_names = c("msp_0005"), gtdb_version = 214.0)
#' get_taxonomy(msp_names = c("msp_0003", "msp_0005"), fields = c("phylum"))
get_taxonomy <- function(msp_names = c(), gtdb_version = 220.0, fields = c()) {
  require(httr)
  require(jsonlite)
  
  if (length(msp_names) == 0) {
    stop("Please provide at least one MSP to request.")
  }
  
  if (!gtdb_version %in% c(214.0, 220.0)) {
    stop("Unsupported GTDB version. Available versions: 214.0, 220.0")
  }
  
  server <- "https://biobanks.gmt.bio/taxo/"
  message(paste("Requesting", server, "..."))
  
  request_code <- httr::status_code(httr::GET(url = server))
  if (request_code == 401) {
    stop("Authentication error. Please check your credentials.")
  }
  
  if (!request_code %in% c(200, 422)) {
    stop(
      sprintf(
        "Connection/authentication error. (Status code: %s). Please check your connection or authorization.",
        request_code
      )
    )
  }
  
  results <- list()
  for (msp in msp_names) {
    request_url <- paste0(server, "?msp=", msp, "&gtb=", gtdb_version, "&format=json")
    
    if (length(fields) > 0) {
      request_url <- paste0(request_url, "&columns=", paste(fields, collapse = ","), ",")
    }
    
    request <- GET(url = request_url)
    code <- status_code(request)
    
    if (code != 200) {
      warning_msg <- switch(
        as.character(code),
        "500" = sprintf("Invalid request for %s. Unknown MSP (Status code: %s)", msp, code),
        "400" = sprintf(
          "Invalid request for %s. Invalid field(s). Use get_fields() or get_fields(msp_name) to list valid fields. (Status code: %s)",
          msp, code
        ),
        sprintf("Invalid request for %s. (Status code: %s)", msp, code)
      )
      warning(warning_msg)
    } else {
      msp_taxonomy <- content(request, "text", encoding = "UTF-8")
      results[[msp]] <- fromJSON(msp_taxonomy, flatten = TRUE)
    }
  }
  
  return(results)
}

#' List available taxonomy fields for an MSP
#'
#' Retrieves and prints the field names that can be queried for a given MSP.
#'
#' @param msp_name A single MSP name (default: \code{"msp_0005"}).
#'
#' @return A character vector of available field names (printed to console).
#' @export
#'
#' @examples
#' get_fields()
#' get_fields("msp_0003")
get_fields <- function(msp_name = "msp_0005") {
  print(names(get_taxonomy(c(msp_name))[[msp_name]]))
}
