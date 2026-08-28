#' cpgdirection: the direction of a CpG's effect on gene expression
#'
#' Assigns a direction to the association between methylation at a CpG site and
#' the expression of a nearby gene, for sites with no measured record, and
#' reports how much each answer is worth.
#'
#' Almost everything is called through the \pkg{data.table}, \pkg{stats},
#' \pkg{utils} and \pkg{tools} namespaces explicitly. What must be imported is
#' the non-standard evaluation symbols, plus the three data.table functions used
#' bare inside \code{[.data.table} expressions, where a \code{::} prefix would
#' be evaluated in the wrong frame. The two annotation packages are imported for their
#' data objects: universal mode reaches them through \code{asNamespace()} rather
#' than \code{::}, which R CMD check cannot see, and declaring them here is both
#' honest and enough to satisfy it.
#'
#' @keywords internal
#' @importFrom data.table := .SD .I .N fifelse fcase setorderv %chin%
#' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
#' @importFrom org.Hs.eg.db org.Hs.eg.db
"_PACKAGE"
