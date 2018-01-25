
knitr::purl(input,output,documentation=2,quiet=T)

knitr::opts_chunk$set(echo = TRUE)#,fig.width=5,fig.height=5,dpi=125)

# knitr::knit_hooks$set(plot = function(x, options) {
#   paste('<figure><img src="',
#         opts_knit$get('base.url'), paste(x, collapse = '.'),
#         '"><figcaption>', options$fig.cap, '</figcaption></figure>',
#         sep = '')
# })


# # A function for generating captions and cross-references
# 
# fig <- local({
#     i <- 0
#     list(
#         cap=function(refName, text, center=FALSE, col="black", inline=FALSE) {
#             i <<- i + 1
#             ref[[refName]] <<- i
#             css_ctr <- ""
#             if (center) css_ctr <- "text-align:center; display:inline-block; width:100%;"
#             cap_txt <- paste0("<span style=\"color:", col, "; ", css_ctr, "\">Figure ", i, ": ", text , "</span>")
#             anchor <- paste0("<a name=\"", refName, "\"></a>")
#             if (inline) {
#                 paste0(anchor, cap_txt)    
#             } else {
#                 list(anchor=anchor, cap_txt=cap_txt)
#             }
#         },
#         
#         ref=function(refName, link=FALSE, checkRef=TRUE) {
#             
#             ## This function puts in a cross reference to a caption. You refer to the
#             ## caption with the refName that was passed to fig$cap() (not the code chunk name).
#             ## The cross reference can be hyperlinked.
#             
#             if (checkRef && !refName %in% names(ref)) stop(paste0("fig$ref() error: ", refName, " not found"))
#             if (link) {
#                 paste0("<A HREF=\"#", refName, "\">Figure ", ref[[refName]], "</A>")
#             } else {
#                 paste0("Figure ", ref[[refName]])
#             }
#         },
#         
#         ref_all=function(){
#             ## For debugging
#             ref
#         })
# })
# ## This chunk replaces the default hook for processing plots. It achieves the purposes,
# ## of laying out auto-numbered captions, but other functionality may be gone.
# 
# library(knitr)
# knit_hooks$set(plot = function(x, options) {
#     sty <- ""
#     if (options$fig.align == 'default') {
#         sty <- ""
#     } else {
#         sty <- paste0(" style=\"text-align:", options$fig.align, ";\"")
#     }
#     
#     if (is.list(options$fig.cap)) {
#         ## options$fig.cap is a list returned by the function fig$cap()
#         str_caption <- options$fig.cap$cap_txt
#         str_anchr <- options$fig.cap$anchor
#     } else {
#         ## options$fig.cap is a character object (hard coded, no anchor)
#         str_caption <- options$fig.cap
#         str_anchr <- ""
#     }
#     
#     paste('<figure', sty, '>', str_anchr, '<img src="',
#         opts_knit$get('base.url'), paste(x, collapse = '.'),
#         '"><figcaption>', str_caption, '</figcaption></figure>',
#         sep = '')
#     
# })
# ## This chucnk will read through *this* Rmd file, and attempt to extract all of the 
# ## labels (not caption text) used for Figure captions. These labels are used
# ## as anchors, so scanning through the document now will allow us to create cross references
# ## before the caption actually appears. 
# 
# ## Get the name of this Rmd file
# rmdFn <- knitr::current_input()  # filename of input document
# 
# ## Read lines and close connection
# rmdCon <- file(rmdFn, open = "r")
# rmdLines <- readLines(rmdCon)
# close(rmdCon)
# 
# ## Pull out all occurences of at least one back tick, followed 
# ## by any number of characters, followed by fig$cap (all on one line)
# figscap_idx <- grep("`+(.*)fig\\$cap", rmdLines)
# rmdLines <- rmdLines[figscap_idx]
# 
# ## Get rid of everything up until the start of the caption label
# ## This presumes the caption label is the first argument of fig$cap()
# ## E.g., fig.cap = fig$cap("my_label", ...)
# rmdLinesSansPre <- sub("(.*)fig\\$cap(.*?)[\"']", "", rmdLines)
# 
# ## Identify everything up until the first quote
# match_data <- regexpr("(.*?)[\"']", rmdLinesSansPre)
# 
# ## Reduce the length by one, because we're not interested in the final quote
# attr(match_data, "match.length") <- attr(match_data, "match.length") - 1
# 
# ## Extract
# fig_labels <- regmatches(rmdLinesSansPre, match_data, invert=FALSE)
# 
# if (length(fig_labels) > 0) {
# 
#     ## Test for duplicates
#     if (anyDuplicated(fig_labels) > 0) stop("Duplicate caption labels detected")
#     
#     ## Create a named list of Figure numbers
#     ref <- as.list(1:length(fig_labels))
#     names(ref) <- fig_labels
# }    
