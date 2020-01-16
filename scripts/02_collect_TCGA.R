### Download GDC data on methylation and gene expression using TCGAbiolinks library.


### Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

if (!require("TCGAbiolinks")) {
  BiocManager::install("TCGAbiolinks")
  }

library(TCGAbiolinks)


### Setup directories and get all projects ids
print(paste("Current directory: ", getwd(), sep=''))
DATA_DIR = "../data/tcga/"
PROJECTS_DIR = paste(DATA_DIR, 'projects/', sep='')

projects = TCGAbiolinks:::getGDCprojects()$project_id


### Define
getProjectData = function(projects) {
  
  # Gets dataframe of samples metadata from queried projects and datatypes.
  # :params
  #   projects --
  #   datatypes --
  # :return
  #   df
  
  dtypes = c("Transcriptome Profiling", "DNA Methylation")
  counts = 0
  i = 0
  
  for (project in projects) {
    
    meta.data = TCGAbiolinks:::getProjectSummary(project)
    
    # If all dtypes not in project, skip
    if (sum(dtypes %in% meta.data$data_categories$data_category)!=2){
      print(project)
      next
    }
    
    
    # Create project folder if doesn't exist
    project_dir = paste(PROJECTS_DIR, project, '/', sep='')
    dir.create(project_dir, showWarnings = FALSE)
    
    
    # Get common patients with both types of data
    query.met = GDCquery(project = project,
                        data.category = "DNA Methylation",
                        legacy = FALSE,
                        platform = c("Illumina Human Methylation 450"))
      
    query.exp <- GDCquery(project = project,
                         # legacy = FALSE,
                          data.category = "Transcriptome Profiling",
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - FPKM-UQ")

    common.patients <- intersect(substr(getResults(query.met, cols = "cases"), 1, 12),
                                 substr(getResults(query.exp, cols = "cases"), 1, 12))
    counts = counts + length(common.patients)
    

    # Query data for the common patients
    query.met <- GDCquery(project = project,
                          #access = 'open',
                          data.category = "DNA Methylation",
                          legacy = FALSE,
                          platform = c("Illumina Human Methylation 450"),
                          barcode = common.patients)
    query.exp <- GDCquery(project = project,
                          #access = 'open',
                          data.category = "Transcriptome Profiling",
                          #legacy = FALSE,
                          data.type = "Gene Expression Quantification", 
                          workflow.type = "HTSeq - FPKM-UQ",
                          barcode = common.patients)
    
    
    # Get and save manifests
    manifest.met = getManifest(query.met, save=FALSE)
    manifest.exp = getManifest(query.exp, save=FALSE)
    
    write.table(manifest.met, paste(project_dir, project, '_manifest_meth.txt', sep=''), 
                quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    write.table(manifest.exp, paste(project_dir, project, '_manifest_expr.txt', sep=''), 
                quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    
    write.table(query.met[,2:11], file=paste(project_dir, project, '_query_meth.tsv',sep=''), 
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(query.exp[,2:11], file=paste(project_dir, project, '_query_expr.csv', sep=''), 
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    write.table(query.met[,1][[1]], file=paste(project_dir, project, '_methylation.csv', sep=''),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    write.table(query.exp[,1][[1]], file=paste(project_dir, project, '_expression.csv', sep=''),
                quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
    
    
    # Download data
    #if ( !dir.exists(paste(project_dir, 'data/methylation/', sep='')) ) {
      #GDCdownload(query.met, method='client', directory=paste(project_dir, 'data/methylation/',sep=''), files.per.chunk=NULL)
    #}
    
    if ( !dir.exists(paste(project_dir, 'data/expression/', sep='')) ) {
      GDCdownload(query.exp, method='client', directory=paste(project_dir, 'data/expression/',sep=''), files.per.chunk=NULL)
    }

    i = i + 1
    if (i==5) {break}
    
  }
  
  # Using GDCDownload downloads gdc-client tool locally, remove it
  file.remove(grep('gdc', list.files(), value=TRUE))
  
  return(counts)
}


### Get data
dat = getProjectData(projects)