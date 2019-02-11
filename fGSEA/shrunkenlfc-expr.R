#!/bigseq/analysis/images/r_3.4.3_v1.simg
library(SummarizedExperiment)


## This script performs differential expression analysis on RNA-seq & Affymetrix data ##
# e.g.
#  ./differential-expr.R --se_obj_file KI_SU-Glom.rnaseq.se-obj.rda --design Group --factor Group --test_name DN --num DN --ref_name HC --den HC
#
#  ./differential-expr.R --se_obj_file rpc2-glomerular.affymetrix.se-obj.rda \
#                        --design disease --factor disease --test_name CKD \
#                        --num Diabetic.Nephropathy,Focal.Segmental.Glomerulosclerosis,Focal.Segmental.Glomerulosclerosis...Minimal.Change.Disease,Hypertension,IgA.Nephropathy,Membranous.Glomerulonephritis,Minimal.Change.Disease,Rapidly.Progressive.Glomerulonephritis,Systemic.Lupus.Erythematous,Thin.Basement.Membrane.Disease,Tumor.Nephrectomy \
#                        --ref_name LD --den Living.Donor


args <- commandArgs( trailingOnly=TRUE )
args

## Read file containing the **serialised** 'SummarizedExperiment' object & check its class is correct ##
if (! "--se_obj_file" %in% args) {
  cat("\n")
  stop("'SummarizedExperiment' flag [--se_obj_file] absent")
} else {
  
  se.file.idx  <- grep("--se_obj_file", args)
  se.file.path <- args[ se.file.idx+1 ]
  
  if (file.exists(se.file.path)) {
    cat("\nLoading _serialised_ 'SummarizedExperiment' object from \"",se.file.path,"\" ... ", sep="")
    summarizedExpt <- readRDS(se.file.path)
    cat("done\n")
    
    check.class <- class(summarizedExpt)[1]
    
    if (check.class=="SummarizedExperiment") {
      cat("\n'SummarizedExperiment' object LOADED :-\n")
      print(summarizedExpt)
      cat("\n")
    }
    else {
      cat("\nObject class = '",check.class,"'\n", sep="")
      stop("Object loaded **IS NOT OF CLASS** 'SummarizedExperiment'")
    }
   
  }
  else {
    cat("\nLocation of file containing 'SummarizedExperiment' object : \"",se.file.path,"\"\n", sep="")
    stop("File **DOES NOT EXIST**")
  }
}
## ------------------------------------------------------------------------- ##


## Check the contrast factor & its associated levels are present ##
## in the 'colData' slot of the 'SummarizedExperiment' object    ##

if ("--test_name" %in% args) { test.id <- args[ grep("--test_name", args)+1 ]
} else { cat("\n"); stop("Name of test levels (Denominator) [--test_name] missing") }

if ("--ref_name" %in% args) { ref.id <- args[ grep("--ref_name", args)+1 ]
} else { cat("\n"); stop("Name of reference levels (Numerator) [--ref_name] missing") }


if ("--factor" %in% args) { contrast.fact <- args[ grep("--factor", args)+1 ]
} else { cat("\n"); stop("Factor name flag [--factor] missing") }

if (contrast.fact %in% colnames(colData(summarizedExpt))) {

  contrast.levs <- as.character( colData(summarizedExpt)[[ contrast.fact ]] )
  cat("\nLevels of contrast factor \"", contrast.fact, "\" : '", paste(sort(unique(contrast.levs)), collapse="', '"),"'\n", sep="")

  if ("--num" %in% args) {
    test.Levels <- unique(unlist( strsplit(args[ grep("--num", args)+1 ], ",", fixed=TRUE) ))
    
    if (! all(test.Levels %in% contrast.levs)) {
      missing.test.levels <- sort( test.Levels[ ! test.Levels %in% contrast.levs ] )
      cat("\nMissing 'Test' [Numerator] levels : '", paste(missing.test.levels, collapse="', '"),"'\n", sep="")
      stop("Levels absent for the contrast factor")
    }
  } else { cat("\n"); stop("Numerator flag [--num] absent") }
  
  if ("--den" %in% args) {
    ref.Levels <- unique(unlist( strsplit(args[ grep("--den", args)+1 ], ",", fixed=TRUE) ))
    
    if (! all(ref.Levels %in% contrast.levs)) {
      missing.ref.levels <- sort( ref.Levels[ ! ref.Levels %in% contrast.levs ] )
      cat("\nMissing 'Reference' [Denominator] level(s) : '", paste(missing.ref.levels, collapse="', '"),"'\n", sep="")
      stop("Level(s) absent for the contrast factor")
    }
  } else { cat("\n"); stop("Denominator flag [--den] absent") }
  
  shared.Levels <- intersect(test.Levels, ref.Levels)
  if (length(shared.Levels)>0) {
    shared.strg <- paste(sort(shared.Levels), collapse="', '")
    cat("\nShared 'Test' & 'Reference' levels : '",shared.strg,"'\n", sep="")
    stop("Common levels")
  }

} else {
  cat("\nContrast factor = '",contrast.fact,"' :-\n", sep="")
  print(colData(summarizedExpt))
  stop("Factor **NOT PRESENT** in sample sheet")
}
## ------------------------------------------------------------------------- ##

# Design formula #
if ("--design" %in% args) { expt.Design <- args[ grep("--design", args)+1 ]
} else { cat("\n"); stop("Design formula flag missing") }



                  ### Differential expression analysis ###

                           ## **RNA-seq** ##

## Function to return the results (logFC & adjusted p-value) of the Wald test ##
waldTest.results <- function(deseq.data.set, factor.name, numerator, denominator, ind.filt=FALSE) {

  IndFilt.flag <- ifelse(ind.filt, "ON", "OFF")
  cat("\nContrast: '",numerator,"' - '",denominator,"' [",factor.name," / Independent Filtering=",IndFilt.flag,"] :-",sep="") 

  wald.test.res <- lfcShrink(deseq.data.set, contrast=c(factor.name, numerator, denominator), independentFiltering=ind.filt)
  summary(wald.test.res)
  cat("\n")
  
  wald.test.vals <- data.frame(
    log2FoldChange=wald.test.res$log2FoldChange, prob=wald.test.res$pvalue,
    row.names=rownames(deseq.data.set)
  )
  return(wald.test.vals)
}


Diff.Expr.NegBin <- function(se.path, se, covariates, design.factor, test.name, num.levels, ref.name, den.levels) {

    if ("dabg" %in% names(rowData(se))) {
      dabg.fail.idx  <- which(! rowData(se)$dabg)
      detected.genes <- rownames(se)[ -dabg.fail.idx ]
      cat("\nOmitting ",nrow(se)-length(detected.genes)," genes with less than 10 counts across *ALL* the samples\n", sep="")
    }
    else { detected.genes <- rownames(se) }
    
    if ("multiQC.pass" %in% colnames(colData(se))) {
      qc.fail.idx <- which(! colData(se)$multiQC.pass)
      if (length(qc.fail.idx)==0) { passed.samples <- rownames(colData(se)) }
      else {
        passed.samples <- rownames(colData(se))[ -qc.fail.idx ]
        cat("\n",length(qc.fail.idx)," samples *FAILED* MultiQC -- omitted from differential expression analysis\n", sep="")
      }
    }
    else { passed.samples <- rownames(colData(se)) }
    
    design.formula <- as.formula( paste0("~",covariates) )
    # See below links for info about 'Reduce'
    # https://stackoverflow.com/questions/28545688/understand-the-reduce-function
    # https://stackoverflow.com/questions/14671172/how-to-convert-r-formula-to-text
    des.form.strg  <- Reduce(paste, deparse(design.formula, width.cutoff=500L))
    des.form.strg  <- gsub("\\s{2,}", " ", des.form.strg, perl=T)
    
    
    cat("\nCreating DESeqDataSet object from 'SummarizedExperiment' object [",length(detected.genes)," \"expressed genes\" / ",length(passed.samples)," samples]\n\n", sep="")
    ddsMat <- DESeqDataSetFromMatrix(
      countData=assays(se)$counts[detected.genes, passed.samples], colData=colData(se)[passed.samples, ], design=design.formula
    )
    print(ddsMat)
    
    cat("\nPerforming differential expression analysis based on the Negative Binomial distribution :-\n")
    dds <- DESeq(ddsMat)
    
    if (length(num.levels)>1) num.strg <- paste0("(",paste(num.levels, collapse=" + "),")/",length(num.levels))
    else                      num.strg <- num.levels
    
    if (length(den.levels)>1) den.strg <- paste0("(",paste(den.levels, collapse=" + "),")/",length(den.levels))
    else                      den.strg <- den.levels
    
    WaldTest.res <- waldTest.results(dds, design.factor, num.strg, den.strg)
    
    de.results.lst <- list(
      se.obj.file=se.path,
      Factor=design.factor, Design.Formula=des.form.strg,
      NumeratorName=test.name, NumeratorLevels=num.levels, DenominatorName=ref.name, DenominatorLevels=den.levels,
      Contrast=paste0(num.strg," - ",den.strg),
      Test.Res=WaldTest.res
    )
    
    dabg.threshold <- 0.5
    cat("\nUsing a detection threshold of <",dabg.threshold," TPM - see \"https://www.ebi.ac.uk/gxa/FAQ.html\"\n\n", sep="")
    
    num.samples <- colnames(se)[ se[[design.factor]] %in% num.levels ]
    num.samples <- num.samples[ num.samples %in% passed.samples ]
    cat("Numerator samples (",length(num.samples),") :-\n", sep=""); print(num.samples)
    
    den.samples <- colnames(se)[ se[[design.factor]] %in% den.levels ]
    den.samples <- den.samples[ den.samples %in% passed.samples ]
    cat("\nDenominator samples (",length(den.samples),") :-\n", sep=""); print(den.samples)
    
    cat("\nFlagging which of the two groups where **ALL** samples are below the threshold .. ")
    dabg.grp.calls <- apply(assays(se)$tpms[ rownames(de.results.lst$Test.Res), ], 1,
      function(row)
      {
        sum( all(row[num.samples]<dabg.threshold), all(row[den.samples]<dabg.threshold) )
      }
    )
    de.results.lst[[ "Test.Res" ]] <- cbind(de.results.lst[["Test.Res"]], fail.detect.cnt=dabg.grp.calls)
    cat("done\n\n")
    
    
    out.file <- paste0(design.factor,"_",test.name,"_vs_",ref.name, ".de-lst.rds")
    cat("Writing differential expression results list ('de.results.lst') to file \"",out.file,"\"  ...  ", sep="")
    saveRDS(de.results.lst, file=out.file)
    cat("finished\n\n")
    return(de.results.lst)
}

## ------------------------------------------------------------------------- ##

                           ## **Affymetrix** ##

Diff.Expr.Affy <- function(se.path, se, covariates, design.factor, test.name, num.levels, ref.name, den.levels) {

  contr.groups   <- colData(se)[[ design.factor ]]
  design.formula <- as.formula( paste0("~0 + ",covariates) )
  
  des.form.strg  <- Reduce(paste, deparse(design.formula, width.cutoff=500L))
  des.form.strg  <- gsub("\\s{2,}", " ", des.form.strg, perl=T)
  
  design <- model.matrix(~0 + contr.groups)
  dimnames(design) <- list(rownames(colData(se)), levels(contr.groups))
  
  if ("dabg" %in% names(rowData(se))) {
    dabg.fail.idx  <- which(! rowData(se)$dabg)
    detected.genes <- rownames(se)[ -dabg.fail.idx ]
    cat("\nOmitting ",nrow(se)-length(detected.genes)," genes : NOT detected above background", sep="")
  }
  else { detected.genes <- rownames(se) }
  
  cat("\nFitting linear model for each retained gene [",length(detected.genes),"]\n\n", sep="")
  rma.norm.fit <- lmFit(assays(se)$normVals[detected.genes,], design)
  
  if (length(num.levels)>1) test.levels <- paste0("(",paste(num.levels, collapse=" + "),")/",length(num.levels))
  else                      test.levels <- num.levels
  
  if (length(den.levels)>1) ref.levels <- paste0("(",paste(den.levels, collapse=" + "),")/",length(den.levels))
  else                      ref.levels <- den.levels
  
  contrast.strg <- paste0(test.levels," - ",ref.levels)
  cat("Contrast matrix :-\n")
  contrasts.matrix <- makeContrasts(contrasts=contrast.strg, levels=design)
  print(contrasts.matrix)
  cat("\n")
  
  cat("Computing coeffs & std err for the contrasts, then compute raw statistics for them ... ")
  rma.contrasts.fit <- eBayes( contrasts.fit(rma.norm.fit, contrasts.matrix) )
  cat("done\n\n")
  
  de.results.lst <- list(
    se.obj.file=se.path,
    Factor=design.factor, Design.Formula=des.form.strg,
    NumeratorName=test.name, NumeratorLevels=num.levels, DenominatorName=ref.name, DenominatorLevels=den.levels,
    Contrast=contrast.strg, 
    Test.Res=setNames( data.frame(
      rma.contrasts.fit$coefficients, rma.contrasts.fit$F.p.value, row.names=rownames(rma.contrasts.fit$coefficients)
    ), c("log2FoldChange", "prob"))
  )
  
  if (exists("dabg.fail.idx"))
  {
    dabg.threshold <- max( assays(se)$normVals[dabg.fail.idx,] )
    cat("\nGenes falling below background in all samples exist!! -- threshold=",dabg.threshold,"\n[Max. value amongst these genes & samples]\n\n", sep="")
    
    num.samples <- colnames(se)[ se[[design.factor]] %in% num.levels ]
    cat("Numerator samples :-\n"); print(num.samples)
    den.samples <- colnames(se)[ se[[design.factor]] %in% den.levels ]
    cat("\nDenominator samples :-\n"); print(den.samples)
    
    cat("\nFlagging which of the two groups where **ALL** samples are below the threshold .. ")
    dabg.grp.calls <- apply(assays(se)$normVals[ -dabg.fail.idx, ], 1,
      function(row)
      {
        sum( all(row[num.samples]<dabg.threshold), all(row[den.samples]<dabg.threshold) )
      }
    )
    de.results.lst[[ "Test.Res" ]] <- cbind(de.results.lst[["Test.Res"]], fail.detect.cnt=dabg.grp.calls)
    cat("done\n\n")
  }

  out.file <- paste0(design.factor,"_",test.name,"_vs_",ref.name, ".shrunken.de-lst.rds")
  cat("Writing differential expression results list ('de.results.lst') to file \"",out.file,"\"  ...  ", sep="")
  saveRDS(de.results.lst, file=out.file)
  cat("finished\n\n")
  return(de.results.lst)
}

## ------------------------------------------------------------------------- ##


if (! "expt.type" %in% names(metadata(summarizedExpt))) {
  cat("\n"); stop( "Experiment type/platform **IS NOT DEFINED** in the 'SummarizedExperiment' object" )
}

if (! metadata(summarizedExpt)$expt.type %in% c("Affymetrix", "RNA-seq")) {
  cat("\n"); stop( "Experiment type **IS NOT** EITHER 'Affymetrix' OR 'RNA-seq'" )
}


if (metadata(summarizedExpt)$expt.type=="RNA-seq") {
  cat("\nLoading BioConductor package 'DESeq2' .. ")
  library(DESeq2)
  cat("done\n")
  
  deTest.res.lst <- Diff.Expr.NegBin(
    se.file.path, summarizedExpt, expt.Design, contrast.fact, test.id, test.Levels, ref.id, ref.Levels
  )
}


if (metadata(summarizedExpt)$expt.type=="Affymetrix") {
  cat("\nLoading BioConductor package 'limma' .. ")
  library(limma)
  cat("done\n")
  
  deTest.res.lst <- Diff.Expr.Affy(
    se.file.path, summarizedExpt, expt.Design, contrast.fact, test.id, test.Levels, ref.id, ref.Levels
  )
}

cat("\n")
str(deTest.res.lst)
cat("\n\n")





















