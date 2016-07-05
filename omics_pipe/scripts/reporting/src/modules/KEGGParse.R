require(KEGGREST)
require(stringr)
require(gdata)
require(AnnotationDbi)

.processENTRY <- function(ENTRY) {
  x <- as.vector(ENTRY)
  return(x)
}

.processDrugNAME <- function(NAME) {
  x <- unlist(strsplit(NAME, split=','))
  country <- sapply(x, function(x) {
    z <- str_match(x, '\\([A-Z]+\\)')
    z <- gsub('\\(', '', z)
    z <- gsub('\\)', '', z)
  })
  drug <- sapply(x, function(x) {
    z <- unlist(strsplit(x, split=' '))
    lenz <- length(z)
    z <- z[-lenz]
    z <- paste(z, collapse=' ')
    return(z)
  })
  names(drug) <- country
}

.processGeneNAME <- function(NAME) {
  x <- unlist(strsplit(NAME, split=','))
  x <- gsub(' ', '', x)
  return(x)
}

.processDRUG_TARGET <- function(DRUG_TARGET) {
  x <- strsplit(DRUG_TARGET, split=':')
  x <- trim(x)
  xx <- lapply(x, function(x) unlist(strsplit(x[2], split=' ')))
  names(xx) <- sapply(x,'[[',1)
  xxx <- unlist2(xx)
  return(xxx)
}

.processCLASS <- function(CLASS) {
  x <- unlist(strsplit(CLASS, split=';'))
  x <- trim(x)
  return(x)
}

.processMOTIF <- function(MOTIF) {
  
}

.processDBLINKS <- function(DBLINKS) {
  x <- strsplit(DBLINKS, split=':')
  xx <- lapply(x, function(x) gsub(' ', '', x))
  xxx <- sapply(xx, "[[", 1)
  names(xxx) <- sapply(xx, "[[", 2)
  return(xxx)
}

.processSTRUCTURE <- function(STRUCTURE) {
  
}

.processACTIVITY <- function(ACTIVITY) {
  x <- strsplit(ACTIVITY, split='\\[')
  lenx <- length(x)
  xx <- lapply(x, function(x) gsub(']', '', x))
  
  ids <- lapply(xx, function(x) {
    z <- unlist(strsplit(x[2], split=' '))
    zz <- substr(z, nchar(z)-5, nchar(z))
    return(zz)
  })
  
  head <- lapply(xx, function(x) {
    z <- x[1]
    return(z)
  })
  head <- unlist(head)
  names(ids) <- head
  return(unlist2(ids))
}

.processTARGET <- function(TARGET) {
  .targetFunc <- function(target) {
    funcDic <- c('inhibitor', 'agonist', 'antagonist', 'blocker', 'activator')
    func <- sapply(target, function(x) {
      res <- sapply(funcDic, function(y) {
        !is.na(pmatch(y, unlist(strsplit(x, split=' '))))
      })
    })
    apply(func, 2, function(x) {
      if(!any(x)) {
        return(NA)
      } else {
        return(names(which(x)))
      }
    })
  }
  if ('TARGET' %in% names(TARGET)) {
    x <- as.vector(sapply(TARGET$TARGET, function(x) str_match(x, 'HSA:[0-9]+')[1]))
    res1 <- gsub('HSA:', '', x)
    targetFunc <- .targetFunc(TARGET$TARGET)
    names(targetFunc) <- res1
  }
  if ('PATHWAY' %in% names(TARGET)) {
    npw <- length(TARGET$PATHWAY)
    pw <- as.vector(sapply(TARGET$PATHWAY, function(x) str_match(x, 'hsa[0-9]+')[1])) # extract the pws
    
    genes <- lapply(TARGET$PATHWAY, function(x) {  # extract the genes for each pw
      z <- str_match(x, '\\([0-9+]+\\)')[1]
      z <- gsub('\\(', '', z)
      z <- gsub('\\)', '', z)
      if (any(grep('+', z))) {
        z <- unlist(strsplit(z, split='\\+'))
      }
      return(z)    
    })
    names(genes) <- pw
    res2 <- unlist2(genes)
  } else {
    x <- as.vector(sapply(TARGET, function(x) str_match(x, 'HSA:[0-9]+')[1]))
    res <- gsub('HSA:', '', x)
    targetFunc <- .targetFunc(TARGET)
    names(targetFunc) <- res
    return(list(TARGET=res, PATHWAY=NA, TARGET_FUNC=targetFunc))
  }
  return(list(TARGET=res1, PATHWAY=res2, TARGET_FUNC=targetFunc)) 
}

.processMETABOLISM <- function(METABOLISM) {
  nmeta <- length(METABOLISM)
  func <- sapply(strsplit(METABOLISM, split=':'), '[[', 1)
  genes <- lapply(METABOLISM, function(x) {
    z <- str_match(x, 'HSA:[0-9+]+')[1]
    z <- gsub('HSA:', '', z)
    if (any(grep('+', z))) {
      z <- unlist(strsplit(z, split='\\+'))
    }
    return(z)    
  })
  names(genes) <- func
  res <- unlist2(genes)
  return(res)
}

.processBRITE <- function(BRITE) {
  ids <- str_match(BRITE, 'BR:br[0-9]+')
  ids <- as.vector(na.omit(ids))
  ids <- gsub('BR:', '', ids)
  return(ids)
}

.processPRODUCT <- function(PRODUCT) {
  setid <- sapply(PRODUCT, function(x) {
    x <- unlist(strsplit(x, split=' '))
    lenx <- length(x)
    x <- x[lenx]
    return(x)
  })
  setid <- as.vector(setid)
  url <- paste('http://dailymed.nlm.nih.gov/dailymed/lookup.cfm?setid=', setid, sep='')
  pname <- sapply(PRODUCT, function(x) {
    x <- unlist(strsplit(x, split=' '))
    x <- x[1]
    return(x)
  })
  pname <- as.vector(pname)
  res <- url
  names(res) <- pname
  return(res)
}

.processSTR_MAP <- function(STR_MAP, KEGGDrugID) {
  map <- sapply(STR_MAP, function(x) {
    x <- unlist(strsplit(x, split=' '))
    x <- x[1]
    return(x)
  })
  url <- paste('http://www.genome.jp/kegg-bin/show_pathway?', map, '+', KEGGDrugID, sep='')
  names(url) <- map
  return(url)
}

.processINTERACTION <- function(INTERACTION) {
  interaction <- sapply(INTERACTION, function(x) {
    x <- unlist(strsplit(x, split=' '))
    lenx <- length(x)
    x <- x[lenx]
    return(x)
  })
  interaction <- as.vector(interaction)
  gene <- sapply(INTERACTION, function(x) {
    x <- str_match(x, 'HSA:[0-9]+')
    x <- gsub('HSA:', '', x)
    return(x)
  })
  names(gene) <- interaction
  return(gene)
}

.processDRUG <- function(DRUG) {
  did <- DRUG[seq(1,length(DRUG), by=2)]
  dname <- DRUG[seq(2,length(DRUG), by=2)]
  dname <- sapply(dname, function(x) {
    z <- unlist(strsplit(x, split=' '))[1]
    return(z)
  })
  names(dname) <- did
  return(dname)
}

.processORGANISM <- function(ORGANISM) {
  z <- str_match(ORGANISM, '(GN:[a-z]+)')[1]
  zz <- gsub('[A-Z]+:', '', z)
  return(zz)
}

.processGENE <- function(GENE) {
  geneid <- GENE[seq(1, length(GENE), by=2)]
  genesymbol <- sapply(GENE[seq(2, length(GENE), by=2)], function(x) {
    z <- unlist(strsplit(x, split=';'))
    zz <- z[1]
    return(zz)
  })
  names(genesymbol) <- geneid
  return(genesymbol)
}

KEGGGetGene <- function(entrezID) {
  tryCatch({
    down <- keggGet(paste('hsa:', entrezID, sep=''))[[1]]
    
    if(is.null(down$ENTRY)) {
      ENTRY <- NA
    } else {
      ENTRY <- .processENTRY(down$ENTRY)
    }  
    if(is.null(down$NAME)) {
      NAME <- NA
    } else {
      NAME <-.processGeneNAME(down$NAME)
    }  
    if(is.null(down$DEFINITION)) {
      DEFINITION <- NA
    } else {
      DEFINITION <- down$DEFINITION
    }  
    if(is.null(down$ORTHOLOGY)) {
      ORTHOLOGY <- NA
    } else {
      ORTHOLOGY <- down$ORTHOLOGY
    }
    if(is.null(down$ORGANISM)) {
      ORGANISM <- NA
    } else {
      ORGANISM <- down$ORGANISM
    }  
    if(is.null(down$PATHWAY)) {
      PATHWAY <- NA
    } else {
      PATHWAY <- down$PATHWAY
    }  
    if(is.null(down$DRUG_TARGET)) {
      DRUG_TARGET <- NA
    } else {
      DRUG_TARGET <- .processDRUG_TARGET(down$DRUG_TARGET)
    }  
    if(is.null(down$POSITION)) {
      POSITION <- NA
    } else {
      POSITION <- down$POSITION
    }    
    if(is.null(down$DBLINKS)) {
      DBLINKS <- NA
    } else {
      DBLINKS <- .processDBLINKS(down$DBLINKS)
    }  
    x <- list(ENTRY=ENTRY,
              NAME=NAME,
              DEFINITION=DEFINITION,
              ORTHOLOGY=ORTHOLOGY,
              ORGANISM=ORGANISM,
              PATHWAY=PATHWAY,
              DRUG_TARGET=DRUG_TARGET,
              CLASS=NA, #.processCLASS(down$CLASS)
              POSITION=POSITION,
              MOTIF=NA, #.processMOTIF(down$MOTIF)
              DBLINKS=DBLINKS,
              STRUCTURE=NA, #.processSTRUCTURE(down$STRUCTURE),
              AASEQ=NA, #ifelse(is.null(down$AASEQ), NA, down$AASEQ),
              NTSEQ=NA #ifelse(is.null(down$NTSEQ), NA, down$NTSEQ)
    )
    return(x)
  }, error = function(e) {
    x <- list(ENTRY=NA,
              NAME=NA,
              DEFINITION=NA,
              ORTHOLOGY=NA,
              ORGANISM=NA,
              PATHWAY=NA,
              DRUG_TARGET=NA,
              CLASS=NA, 
              POSITION=NA,
              MOTIF=NA, 
              DBLINKS=NA,
              STRUCTURE=NA,
              AASEQ=NA, 
              NTSEQ=NA 
              )
    return(x)
  }
  )
}


KEGGGetDrug <- function(KEGGDrugID) {
  tryCatch({
    down <- keggGet(KEGGDrugID)[[1]]
    
    if(is.null(down$ENTRY)) {
      ENTRY <- NA
    } else {
      ENTRY <- .processENTRY(down$ENTRY)
    }
    if(is.null(down$NAME)) {
      NAME <- NA
    } else {
      NAME <- .processDrugNAME(down$NAME)
    } 
    if(is.null(down$FORMULA)) {
      FORMULA <- NA
    } else {
      FORMULA <- down$FORMULA
    } 
    if(is.null(down$EXACT_MASS)) {
      EXACT_MASS <- NA
    } else {
      EXACT_MASS <- down$EXACT_MASS
    }
    if(is.null(down$MOL_WEIGHT)) {
      MOL_WEIGHT <- NA
    } else {
      MOL_WEIGHT <- down$MOL_WEIGHT
    }  
    if(is.null(down$MOL_WEIGHT)) {
      MOL_WEIGHT <- NA
    } else {
      MOL_WEIGHT <- down$MOL_WEIGHT
    }  
    if(is.null(down$ACTIVITY)) {
      ACTIVITY <- NA
    } else {
      ACTIVITY <- .processACTIVITY(down$ACTIVITY)
    }
    if(is.null(down$REMARK)) {
      REMARK <- NA
    } else {
      REMARK <- down$REMARK
    }    
    if(is.null(down$COMMENT)) {
      COMMENT <- NA
    } else {
      COMMENT <- down$COMMENT
    }    
    if(is.null(down$TARGET)) {
      TARGET <- list(TARGET=NA, PATHWAY=NA, TARGET_FUNC=NA)
    } else {
      TARGET <- .processTARGET(down$TARGET)
    }    
    if(is.null(down$METABOLISM)) {
      METABOLISM <- NA
    } else {
      METABOLISM <- .processMETABOLISM(down$METABOLISM)
    } 
    if(is.null(down$INTERACTION)) {
      INTERACTION <- NA
    } else {
      INTERACTION <- .processINTERACTION(down$INTERACTION)
    }  
    if(is.null(down$STR_MAP)) {
      STR_MAP <- NA
    } else {
      STR_MAP <- .processSTR_MAP(down$STR_MAP, KEGGDrugID)
    }  
    if(is.null(down$PRODUCT)) {
      PRODUCT <- NA
    } else {
      PRODUCT <- .processPRODUCT(down$PRODUCT)
    }  
    if(is.null(down$BRITE)) {
      BRITE <- NA
    } else {
      BRITE <- .processBRITE(down$BRITE)
    }  
    if(is.null(down$DBLINKS)) {
      DBLINKS <- NA
    } else {
      DBLINKS <- .processDBLINKS(down$DBLINKS)
    }  
    if(is.null(down$ATOM)) {
      ATOM <- NA
    } else {
      ATOM <- down$ATOM
    }  
    if(is.null(down$BOND)) {
      BOND <- NA
    } else {
      BOND <- down$BOND
    }  
    x <- list(ENTRY=ENTRY,
              NAME=NAME,
              FORMULA=FORMULA,
              EXACT_MASS=EXACT_MASS,
              MOL_WEIGHT=MOL_WEIGHT,
              ACTIVITY=ACTIVITY,
              REMARK=REMARK,
              COMMENT=COMMENT,
              TARGET=TARGET,
              METABOLISM=METABOLISM,
              INTERACTION=INTERACTION,
              STR_MAP=STR_MAP,
              PRODUCT=PRODUCT,
              BRITE=BRITE,
              DBLINKS=DBLINKS,
              ATOM=ATOM,
              BOND=BOND
              )
    return(x)  
  }, error = function(e) {
    x <- list(ENTRY=NA,
              NAME=NA,
              FORMULA=NA,
              EXACT_MASS=NA,
              MOL_WEIGHT=NA,
              ACTIVITY=NA,
              REMARK=NA,
              COMMENT=NA,
              TARGET=NA,
              METABOLISM=NA,
              INTERACTION=NA,
              STR_MAP=NA,
              PRODUCT=NA,
              BRITE=NA,
              DBLINKS=NA,
              ATOM=NA,
              BOND=NA
              )
    return(x)  
  }
 )
}

KEGGGetPathway <- function(KEGGPwID) {
  tryCatch({
    down <- keggGet(KEGGPwID)[[1]]
    if(is.null(down$ENTRY)) {
      ENTRY <- NA
    } else {
      ENTRY <- .processENTRY(down$ENTRY)
    } 
    if(is.null(down$NAME)) {
      NAME <- NA
    } else {
      NAME <- down$NAME
    } 
    if(is.null(down$DESCRIPTION)) {
      DESCRIPTION <- NA
    } else {
      DESCRIPTION <- down$DESCRIPTION
    } 
    if(is.null(down$CLASS)) {
      CLASS <- NA
    } else {
      CLASS <- .processCLASS(down$CLASS)
    } 
    if(is.null(down$PATHWAY_MAP)) {
      PATHWAY_MAP <- NA
    } else {
      PATHWAY_MAP <- down$PATHWAY_MAP
    } 
    if(is.null(down$DISEASE)) {
      DISEASE <- NA
    } else {
      DISEASE <- down$DISEASE
    } 
    if(is.null(down$DRUG)) {
      DRUG <- NA
    } else {
      DRUG <- .processDRUG(down$DRUG)
    } 
    if(is.null(down$ORGANISM)) {
      ORGANISM <- NA
    } else {
      ORGANISM <- .processORGANISM(down$ORGANISM)
    } 
    if(is.null(down$GENE)) {
      GENE <- NA
    } else {
      GENE <- .processGENE(down$GENE)
    } 
    if(is.null(down$COMPOUND)) {
      COMPOUND <- NA
    } else {
      COMPOUND <- down$COMPOUND
    } 
    if(is.null(down$KO_PATHWAY)) {
      KO_PATHWAY <- NA
    } else {
      KO_PATHWAY <- down$KO_PATHWAY
    } 
    if(is.null(down$REFERENCE)) {
      REFERENCE <- NA
    } else {
      REFERENCE <- down$REFERENCE
    } 
    x <- list(ENTRY=ENTRY,
              NAME=NAME,
              DESCRIPTION=DESCRIPTION,
              CLASS=CLASS,
              PATHWAY_MAP=PATHWAY_MAP,
              DISEASE=DISEASE,
              DRUG=DRUG,
              ORGANISM=ORGANISM,
              GENE=GENE,
              COMPOUND=COMPOUND,
              KO_PATHWAY=KO_PATHWAY,
              REFERENCE=REFERENCE
              )
    return(x)
  }, error = function(e) {
    x <- list(ENTRY=NA,
              NAME=NA,
              DESCRIPTION=NA,
              CLASS=NA,
              PATHWAY_MAP=NA,
              DISEASE=NA,
              DRUG=NA,
              ORGANISM=NA,
              GENE=NA,
              COMPOUND=NA,
              KO_PATHWAY=NA,
              REFERENCE=NA
    )
    return(x)
  }
  )
}




