naCumInc <- function(data,...) 
{
    ## 'data' is either a data.frame containing the variables indicated by 
    ##    'futime','eventVar','groupVar' and 'strataVars'
    ##
    ##   _OR_ 
    ##
    ## it's a object of class 'naCumInc' - the output of an earlier call to this
    ## function), which will then be marginalized over stratification variables
    ##
    ## NOTE:  Only data.frame method is current implemented
    UseMethod("naCumInc", data)
}


## Function to censor data to a particular time-point in follow-up time
censorData <- function(data, censorTime, futimeVar, eventVar, 
                       censInd = FALSE, censIndName="censored") {

  ## 'data' is a data.frame containing the variables whose names are given by 
  ##    variables 'futimeVar' and 'eventVar'
  ##
  ## 'censorTime is the time to censor the data to 
  ##
  ## 'futimeVar' is a character string giving the name of the variable
  ##     containing the follow-up time 
  ##
  ## 'eventVar'  is a character string giving the name of the variable
  ##     containing the event indicator
  ##
  ## 'censInd' is a logical indicating whether the data returned should include
  ##     an indicator of which records were censored
  ##
  ## 'censIndName' is a character string specifying the name to give the censoring
  ##     indicator variable

  ## 'eps' = epsilon = small number used to avoid issues in comparisons of 
  ##  floating-point numbers
  eps <- sqrt( .Machine$double.eps )

  ## Find data points that need to be censored
  censIndVec <-  data[[ futimeVar ]] > (censorTime + eps) 

  ## all censored records have futime = censorTime and event = 0 (if eventVar
  ## is not a factor) or event = first-factor-level
  data[[ futimeVar ]][ censIndVec ] <- censorTime
  newValue <- ifelse( !is.factor(data[[eventVar]]), 0, levels(data[[eventVar]])[1] )
  data[[ eventVar ]][ censIndVec ]  <- newValue

  if ( censInd ) 
     data[[ censIndName ]] <- as.integer(censIndVec)

  return( data )
}


# Future expansion possibility, not currently implemented
#
#naCumInc.naCumInc <- function(data, strataWeights=NULL) {
## code to deal with marginalizing
## If inherits(data, 'naCumInc') then <what?>  
##    e.g. "Then you must be combining over levels of strata and 'strataWeights' must be specified'
##    (not sure yet that this is true)
#} 

## naCumInc.data.frame is the naCumInc method for data.frames (currently the only method)
## naCumInc computes cumulative hazards and cumulative incidence curves based on the
## Nelson-Aalen ('na') estimator of the cumulative hazard. 
naCumInc.data.frame <- 
  function(data,  futimeVar, eventVar,  groupVar=NULL,  strataVar= NULL, idVar=NULL,
           censor=list(minAtRisk= NULL, time= NULL), confLevel=0.95,
           stratifiedEstimator= FALSE, strataWeights= NULL )
{

  # futimeVar: a character string giving the name of a variable in 'data' containing
  #     the follow-up times
  #
  # eventVar: a character string giving the name of a variable in 'data' containing
  #     an integer vector indicating the occurence of an event (1=yes, 0=no)
  #
  # groupVar: a character string giving the name of a variable in 'data'. 
  #    The variable must be categorical and should identify some characteristic of interest.
  #    Typically this would be something like treatment group, but could be anything.
  #
  # strataVar:  a character *string* (i.e. must be length 1), giving the name of the
  #    variable in 'data' to be stratified on.  If you have multiple variables you want to
  #    stratify on (e.g. sex and dosage) you can do this by first creating a single variable
  #    in your dataset that jointly specifies the levels of all those variables and then
  #    pass that variable name via 'strataVar'.  You can create such a variable using 
  #    e.g. the R function interaction(). 
  #
  # idVar: a character string providing the name of the variable in 'data' containing a
  #    unique record identifier.  Can be useful operationally, but not otherwise required 
  #
  # confLevel: gives the desired confidence level for the confidence intervals for: 
  #    cumhaz, cuminc, and cir/PE if requested 
  # 
  # stratifiedEstimator:  a logical (TRUE/FALSE) indicating whether a 'stratified estimator'
  #    should be computed and returned.  If TRUE, then a 'strataVar' must have been specified
  #    and 'strataWeights' provided
  #
  # strataWeights: a NAMED vector of weights to use in computing a stratified estimator
  #    The names must match the levels of the strata variable given by 'strataVar'.
  #    In the case that each group has the same strata levels (i.e. the typical case), then
  #    the weights must sum to 1.  If instead there are sets of strata values applying to
  #    different groups, then the weights must sum to 1 for each set of strata values.   
  #
  # censor: allows a way to specify if/how censoring should be done.  If 'censor' is 
  #    specified, it must be a list with components named either 'minAtRisk' or 'time'.
  #    If 'minAtRisk' = m is specified, then that instructs that censoring be done to 
  #    ensure that the number of participants at risk at all timepoints is >= m.  This is
  #    implemented by finding the latest timepoint such that all group x strata subcohorts
  #    satisfy the 'minimum-at-risk' criteria at that timepoint - that timepoint will be 
  #    returned as 'censorTime'.  If 'time' = t is specified, then the data will be censored
  #    at time 't'.  If *both* 'minAtRisk' and 'time'=t are specified, the time that satisfies
  #    the 'minAtRisk' criteria will be determined, and then the *smaller* of that time and
  #    time 't' will be used for censoring (and *that* value returned via 'censorTime')


  ## Plan:
  ##   Since suvfit is going to take our group and strata variables and make one big
  ##   interaction out of them - which we'll need to decode - it's better if *WE* 
  ##   create the interaction in our data.  Doing this means we can easily link to pick
  ##   up the values of the original variables, which we'll need to do.


  ##  'eps' (epsilon) = a small number used to avoid issues of floating point arithmetic
  eps <- sqrt( .Machine$double.eps )


  ## check all "Var" inputs - must be character vectors and be in 'data'
  if ( !is.null(groupVar) && ( length(groupVar) > 1 || !is.character(groupVar) ) )
    stop("Error: argument 'groupVar' must be a character vector of length 1\n\n")

  if ( !is.null(strataVar) && ( length(strataVar) > 1 || !is.character(strataVar) ) )
    stop("Error: argument 'strataVar' must be a character vector of length 1\n\n")

  if ( !is.null(idVar) && ( length(idVar) > 1 || !is.character(idVar) ) )
    stop("Error: argument 'idVar' must be a character vector of length 1\n\n")

  if ( length(futimeVar) > 1 || !is.character(futimeVar) )
    stop("Error: argument 'futimeVar' must be a character vector of length 1\n\n")

  if ( length(eventVar) > 1 || !is.character(eventVar) )
    stop("Error: argument 'eventVar' must be a character vector of length 1\n\n")


  if (!is.null(groupVar) && ! groupVar %in% names(data) )
    stop("The variable ", groupVar, "was not found in dataset ", data, "\n\n")

  if (! futimeVar %in% names(data) )
    stop("The variable ", futimeVar, "was not found in dataset ", data, "\n\n")

  if (! eventVar %in% names(data) )
    stop("The variable ", eventVar, "was not found in dataset ", data, "\n\n")

  if (!is.null(strataVar) && !strataVar %in% names(data) )
    stop("The variable ", strataVar, "was not found in dataset ", data, "\n\n")

  if (!is.null(idVar) && !idVar %in% names(data) )
    stop("The variable ", idVar, "was not found in dataset ", data, "\n\n")


  if ( !missing(censor) ) {
    if ( !is.list(censor) )
      stop("Error: argument 'censor' must be a list\n\n")

    censor.components <- c('minAtRisk', 'time')
    if ( !all( names(censor)  %in% censor.components ) )
      stop("Error: 'censor' may only have compents with names among: ",
           paste(censor.components, collapse=", "),"\n\n")
  }


  if ( stratifiedEstimator && is.null(strataVar) )
    stop("A stratified estimator was requested, but no value of 'strataVar' provided\n\n")

  if ( stratifiedEstimator && is.null(strataWeights) )
    stop("A stratified estimator was requested, but no value of 'strataWeights' provided\n\n")

  if ( !stratifiedEstimator && !is.null(strataWeights) )
    cat("\n\n", sep="",
        "NOTE: argument strataWeights was specified, but will be ignored.\n",
        "The argument is not relevant unless stratifiedEstimator=TRUE, which it is not\n\n")



  ## Now that we're done checking inputs, let's define our output object structure and 
  ## fill in some pieces
  naObj <- list(
              callArgs = list(
                  futimeVar = futimeVar,
                  eventVar  = eventVar,
                  groupVar  = groupVar,
                  strataVar = strataVar,
                  idVar = idVar,
                  stratifiedEstimator = stratifiedEstimator,
                  strataWeights = strataWeights,
                  censor = censor,
                  confLevel = confLevel
              )[ c(TRUE, TRUE, !is.null(groupVar), !is.null(strataVar), !is.null(idVar),
                   TRUE, stratifiedEstimator, !missing(censor), TRUE ) ],


              ## indicator of whether any censoring was *PERFORMED* on the data
              ## (as opposed to having simply specifyied censoring arguments - since
              ##  it's possible that they'd have no effect)
              censored = FALSE,

              ## summary info on censoring - will be filled in if censoring was ATTEMPTED
              censorInfo = list(
                  censorTime = NA,
                  nRecsCensored = NA,
                  nEventsCensored = NA
                ), 

              ## the data used in estimation is returned here, will be censored if 
              ## output component 'censored' is TRUE
              data = NULL,

              ## data.frame of cumulative incidence and cumulative hazard estimates,
              ## standard errors and confidence intervals, possibly for separate strata, 
              ## e.g. if 'strataVar' was specified but a stratified estimator was not
              ## requested 
              cuminc = NULL 
            )

  ## put a class on the object (the "na" = "nelson-aalen")
  class(naObj) <- c("naCumInc", class(naObj)) 


  ## post-process censor
  if ( missing(censor) ) {

    ## if censor wasn't specified, set it to NULL
    censor <- NULL

  } else {
    ## if it was specified, remove NULL components (which could be both of them!)

    null.components <- sapply(censor, is.null) 

    ## if they're all NULL, set censor to NULL
    if ( all( null.components ) ) {
       censor <- NULL
    } else {
      ## Keep  non-NULL components
      censor <- censor[ !null.components ] 
    }
  }

  ## subset the data to the variables we need, prior to modifying
  d <- data[, c(idVar, futimeVar, eventVar, groupVar, strataVar) ]


  if (is.null(groupVar) && is.null(strataVar) ) {

    ## 'Strata' is created for internal use
    Strata <- FALSE
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ 1")

  } else {
    ## create our own temporary interaction variable, combining 'groupVar' and
    ## 'strataVar'.  This will allow us to unwind some of the changes that will
    ## be made to the levels in the output of survfit.
    ## The code below works whether or not 'strataVar' is specified
    tmpStrataVar <- "my_tmp_Strata"
    d[[tmpStrataVar]] <- interaction(d[,c(groupVar,strataVar)], lex.order=TRUE)

    ## store the mapping of the new strata var. to its component ones
    strataMapping <- unique(d[, c(tmpStrataVar, groupVar, strataVar)] )

    Strata <- TRUE
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ ", tmpStrataVar)
  }


  ## Now we're ready to do censoring, if any was requested
  if ( !is.null(censor) ) { 

    if ( is.null( censor$minAtRisk ) ) {
      ## if minAtRisk not given, then we have only censor$time to work with - EZ stuff!
      censorTime <- censor$time

    } else {     
      ## determine the censorTime corresponding to minAtRisk 

      ## coerce minAtRisk to an integer - in case user was trying to be cute
      ## And make sure it's not a negative number
      minAR <- max(0L, ceiling( censor$minAtRisk ) )

      ## use 'minAR' to determine the time to censor at
      if (!Strata) {
        # if no groups/strata, use all data 
        censorTime <- sort( d[[ futimeVar ]], decreasing=TRUE )[ minAR ] 

      } else {

        # else do by strata
        censTimes.strata <- 
            by(data = d, INDICES = d[[tmpStrataVar]], 
               FUN = function(x,n) sort(x[[ futimeVar ]], decreasing=TRUE )[n],
               n = minAR )

        censorTime <- min( as.vector(censTimes.strata), na.rm=TRUE )
      }

      ## if censor$time was *also* specified, then the minimum of it and the time just computed
      if ( !is.null( censor$time ) ) 
        censorTime <- min( censorTime, censor$time, na.rm=TRUE )
    }

    ## Now do the censoring the data at 'censorTime'
    d <- censorData( 
            data = d,
            censorTime= censorTime, 
            futimeVar= futimeVar, 
            eventVar= eventVar,
            censInd = TRUE ) 

    ## set censoring indicator in output object
    naObj$censored <- ifelse( sum(d$censored) > 0, TRUE, FALSE ) 

    ## set the censoring time that was used (regardless of whether it had an effect)
    naObj$censorInfo$censorTime <- censorTime

    ## if censoring was done, modify data, compute some info and store
    if ( naObj$censored ) {

      ## also attach the original data futime and event indicator columns to the
      ## censored dataset, storing them in variables named by adding a "orig_" prefix
      ## to the original variable name
      d[[ paste0("orig_", futimeVar) ]] <- data[[ futimeVar ]]
      d[[ paste0("orig_", eventVar) ]]  <- data[[ eventVar ]]

      naObj$censorInfo$nRecsCensored   <- sum(d$censored)
      naObj$censorInfo$nEventsCensored <- sum( d[[ paste0("orig_", eventVar) ]][ d$censored == 1 ] )

    } else {

      ## no censoring done
      naObj$censorInfo$nRecsCensored   <- 0
      naObj$censorInfo$nEventsCensored <- 0

      ## remove censoring indicator from data
      d$censored <- NULL
    } 
  }
  ## store data into output object (excluding the temporary strata var we created, if present)
  naObj$data <- d[ , names(d) != tmpStrataVar ] 



  ## load survival package
  library(survival)

  ## coerce our character string 'formula.char' to an actual "formula"
  SurvFormula <- as.formula( formula.char )

  ## Now comes the magic:
  ##
  ## 1. using type="fleming-harrington" and error= "tsiatis" in survift will provide
  ##    us with a survival curve based on the nelson-aalen estimator (i.e. it'll give
  ##    us values of 'surv' computed as: exp(-nelsonAalen) )
  ##
  ## 2. calling summary() on the surfit output with option 'extend=TRUE' and with 
  ##    the 'times' argument given *all* event times (from all groups and strata)
  ##    will provide us with estimates evaluated at the same set of timepoints for 
  ##    every group and strata level.  That means we don't need to worry about that
  ##    issue later when we compute stratified estimates or ratios.
  all.uniq.event.times <- c(0, sort( unique( d[[ futimeVar ]][ d[[eventVar]] == 1 ] ) ) )

  summ <- summary( survfit( SurvFormula, data=d, 
                            type= "fleming-harrington", error= "tsiatis", conf.int=confLevel),
                   extend = TRUE, times = all.uniq.event.times )

  summ.cols <- c("time", "n.risk", "n.event", if (Strata) "strata", 
                 "surv", "std.err", "lower", "upper")
  df.summ   <- as.data.frame( summ[ summ.cols ] )


  ## transform to get cuminc, cumhaz and se.cumhaz, as well as confidence intervals
  ## from 'surv','std.err','upper' and 'lower' ('lower'/'upper' are confint for 'surv')
  df.summ <- within( df.summ, {

               cuminc    <- 1 - surv
               se.cuminc <- std.err
               cumhaz    <- -log(surv)
               se.cumhaz <- std.err/surv
               lo.cumhaz <- -log(upper)
               up.cumhaz <- -log(lower)
               lo.cuminc <- 1 - upper
               up.cuminc <- 1 - lower
             })


  ## drop columns no longer needed: surv, std.err, upper, lower 
  dropCols <- c("surv","std.err","lower","upper")
  cumInc <- df.summ[ , ! names(df.summ) %in% dropCols ]

  ## Remaining columns are: 
  ##  "time", "n.risk", "n.event", maybe "strata", 
  ##  "cuminc", "se.cuminc", "cumhaz", "se.cumhaz"
  ##  "lo.cumhaz", "up.cumhaz", "lo.cuminc", "up.cuminc"

  if ( Strata ) {
    ## fix the levels of 'strata':
    ## survfit alters the level names of 'strata', our code (below) changes them back
    levels( cumInc$strata ) <- sub( paste0(tmpStrataVar,"="), "", levels(cumInc$strata) )

    ## rename 'strata' before merging to add original variables back on (in case one
    ## of the original variables has that same name)
    names(cumInc)[ names(cumInc) == "strata" ] <- tmpStrataVar

    ## deconvolute the 'Strata' variable to extract 'group, and then drop the temp variable
    cumInc <- merge( strataMapping, cumInc ) #, by= get("tmpStrataVar") )

    ## resort the data, by tmpStrataVar and time
    cumInc <- cumInc[ order(cumInc[[tmpStrataVar]], cumInc$time), ] 

    ## drop the tmp strata variable
    cumInc[[ tmpStrataVar ]] <- NULL
  }

  ## re-org the columns (and exclude ones we don't want to keep)
  ordCols <- c(groupVar, if (Strata) strataVar , "time", "n.risk","n.event",
               "cuminc", "se.cuminc", "lo.cuminc", "up.cuminc",
               "cumhaz", "se.cumhaz", "lo.cumhaz", "up.cumhaz" )
  #print( ordCols[ ! ordCols %in% names(df.summ) ] )
  cumInc <- cumInc[ , ordCols ]

  ## if a stratified estimator was NOT requested, then we're done here
  if (!stratifiedEstimator) {
     ## store estimates in output object
     naObj$cuminc <- cumInc

     return( naObj )
  }

  ## if a stratified estimator *was* requested, then we continue on...

  ## get vector of groups to loop over
  if ( is.factor(cumInc[[ groupVar ]] ) ) {
    grps <- levels(cumInc[[ groupVar ]])
  } else {
    grps <- sort(unique(cumInc[[ groupVar ]]))
  }

  ## create a list to store new data into, during the loop
  tmpList <- vector("list", length(grps) )


  ## loop over the levels of groupVar 
  for ( g in grps ) {

    ## identify the rows that correspond to the current group 'g' and extract them
    dat.g <- cumInc[ which( cumInc[[ groupVar ]] == g ), ]

    ## get vector of the strata present in group 'g' 
    strata  <- unique( dat.g[[ strataVar ]] )


    ## if there's more than one strata level, then begin creation of the 
    ## stratified estimator for this group
    if ( length(strata) > 1) {

      ## check if all strata have weights in strataWeights vector
      strataMissingWeights <- strata[ ! strata %in% names(strataWeights) ]

      if ( length( strataMissingWeights ) > 0 )
         stop("strataWeights values were not found for the following strata: \n",
              paste( strataMissingWeights, collapse= " "), "\n",
              "Please ensure that you have named the components of 'strataWeights' with names\n",
              "identical to those in the data (shown here)\n\n")

      ## pull out the strataWeights that correspond to our current strata   
      ## (we're indexing the weights by their names - or trying to...)
      wts <- strataWeights[ strata ]

      ## if weights don't sum to 1 then issue error and bail out
      if ( abs( sum(wts) - 1 ) > eps ) 
        stop("Error in creating stratified estimator for group ", g, ". \n\n",
             "The 'strataWeights' corresponding to the strata levels found for this group,\n",
             "which are: ", paste(strata, collapse= " "), "\n",
             "do not sum to 1 as is required.\n",
             "If the strata levels listed here do not include all that are expected, then this\n",
             "group's data are lacking this level\n\n" )
             
      ## split the data by strata levels, with each level going into a 
      ## diff. component of a list 
      dList <- split( dat.g[, names(dat.g) != strataVar ],  
                      dat.g[[ strataVar ]] )

      ## pull out one element of dList to use to store our stratified estimates into
      ## (we'll overwrite the current variables, as needed)
      stratEst <- dList[[1]]

      ## the 'n.risk' and 'n.event' columns become meaningless for stratified estimators
      ## (or so it would seem), so set these values to NA.  
      ## We keep these columns for consistency of output structure with non-stratified analyses
      stratEst$n.risk  <- NA
      stratEst$n.event <- NA


      ## Now compute stratified 'cumhaz' estimate:
      # --------------
      # 1. extract 'cumhaz' from the data subsets and multiply by the weight
      weightedCumHazList <- lapply( 1:length(dList),
                                    function(i, dL, w) w[i]*dL[[i]]$cumhaz ,
                                    dL = dList, w = wts )

      # 2. cbind the weighted cumhaz together and then sum the rows, store into
      #    'cumhaz' element of stratEst
      stratEst$cumhaz <- rowSums( do.call(cbind, weightedCumHazList) )


      ## Compute the variance of the stratified 'cumhaz' estimate:
      # --------------
      # 1. extract 'se.cumhaz' from the data subsets, multiply by the weight,
      #    and square that
      weightedVarList <- lapply( 1:length(dList),
                                 function(i, dL, w) (w[i]*dL[[i]]$se.cumhaz)^2 ,
                                 dL = dList, w = wts )

      # 2. cbind the weighted var.s together, sum the rows, and then take the
      #    sqrt() to get new 'se.cumhaz'
      stratEst$se.cumhaz <- sqrt( rowSums( do.call(cbind, weightedVarList) ) )


      ## Now compute 'cuminc' and 'se.cuminc' from 'cumhaz' and 'se.cumhaz'
      #stratEst <- transform( stratEst, 
      #                       cuminc    = 1 - exp(-cumhaz),
      #                       se.cuminc = exp(-cumhaz)*se.cumhaz )

      stratEst <- within( stratEst, {
                    cuminc    <- 1 - exp(-cumhaz)
                    se.cuminc <- exp(-cumhaz)*se.cumhaz 
                    lo.cumhaz <- pmax(cumhaz - qnorm(0.975)*se.cumhaz, 0)
                    up.cumhaz <- pmin(cumhaz + qnorm(0.975)*se.cumhaz, 1)
                    lo.cuminc <- 1 - exp(-lo.cumhaz) 
                    up.cuminc <- 1 - exp(-up.cumhaz) 
                  })

      ## Store the estimates
      tmpList[[ which( grps == g) ]] <- stratEst

    } else {
      if ( length( strata ) == 1 ) {

        ## It's possible that there should have been more than one stratum, do a check 
        ## then issue a message  - or an error, to the user

        ## see if a strataWeight value was provided for the strata, and if so, if it's '1'
        wt <- strataWeights[ names(strataWeights) == strata ]
        if ( length(wt) == 1 ) {
          ## if wt isn't equal to 1, issue a warning and quit, otherwise continue on quietly
          if ( abs( wt - 1 ) > eps ) 
            stop(
              "Error in creating stratified estimator for group ", g, ". \n\n",
              "Only one stratum was found for this group, which was ", strata, "\n",
              "and the weight provided for it in 'strataWeights' was not equal to 1\n",
              "Perhaps you were expecting this group to contain more strata?  Exiting.\n\n" )
        } else {

          ## no weight was found, issue informational message then continue on
          cat( "\n\n*************************************************************************\n",
               "*** Important Note ***\n",
               "You have requested that stratified estimators be computed, but only one\n",
               "strata has been found for group ", g, ", which is: ", strata, "\n",
               "Please ensure that this is as expected. Continuing on... \n",
               "**************************************************************************\n\n",
               sep="" )
        }

        ## store the data - minus the strata variable
        tmpList[[ which( grps == g) ]] <- dat.g[ , names(dat.g) != strataVar ]
      }
    }
  }  # end loop over 'grps'
  

  ## combine data from the separate groups together again, and store to output object
  naObj$cuminc <- do.call(rbind, tmpList )

  ## .... and we're done!
  return( naObj )
}



# --------------------------------------------------------------------------------
# - NOTE: version 1 
#   only works if all groups present in the same input object 
# --------------------------------------------------------------------------------

## Efficacy estimation via Cumulative Incidence Ratio (CIR) 

EffCIR <- function(obj, refLvl, cmpLvl=NULL, confLevel=0.95, nullHypEff=NULL) {

  # obj : an object of class 'naCumInc' containing estimates
  #       for the specified values of 'refLvl' and (if given) 'cmpLvl'.
  #
  # refLvl : the "level" of the group variable to be used as the reference 
  #         (i.e. denominator) when computing the cumulative incidence ratio(s).
  #         This variable must be specified.
  #
  # cmpLvl : the level of the group variable to compare to the reference.  Used
  #          as in numerator when constructing the cumulative incidence ratio(s).
  #          If cmpLvl is left as 'NULL', then CIRs will be computed for *all*
  #          non-reference-group levels vs. the specified reference level.
  #
  # confLevel: the confidence level used to construct pointwise intervals for the
  #            CIR and Efficacy estimates
  #
  # nullHypEff:  The value of the null hypothesis that the user wishes to test.  It must
  #    be a numeric vector with all values < 1, more than one hypothesis can be specified.
  #    The values provided should specify the null for the *efficacy* parameter (i.e.
  #    1 - Cumulative-Incidence-Ratio), not for e.g. log(CIR).
  #    If the argument is not specified by the user, the function returns results for
  #    nullHypEff=0 (Efficacy = 0, CIR = 1).
  #
  #    The test provided is the two-sided Wald test of
  #       H0: Efficacy=xxx vs. HA: Efficacy != xxx, 
  #    it is NOT testing the directional hypotheses 
  #       H0: Efficacy <= xxx  vs. HA: Efficacy > xxx


  if ( missing(obj) || !inherits(obj, c("naCumInc","crCumInc")) ) 
    stop("You must provide an object of type 'naCumInc' or 'crCumInc' as the ",
         "first argument to this function.\n\n")

  if ( missing(refLvl) )
    stop("You must specify the argument 'refLvl' \n\n")

  if ( length(refLvl) > 1 || !is.character(refLvl) )
    stop("Argument 'refLvl' must be a character string (a char. vector of length 1) \n\n")

  if ( is.null( obj$callArgs$groupVar ) )
    stop("The object you provided does not contain a group-variable.\n",
         "Please provide an object containing cumulative-incidence curves\n",
         "created by specifying a value for argument 'groupVar'\n")

  ## Do basic check on nullHypEff if it was specified
  if ( !is.null(nullHypEff) ) {
    if (!is.numeric(nullHypEff) || any(nullHypEff >= 1) )
      stop("Argument 'nullHyp' must be a numeric vector with all values < 1 \n\n")
  }


  ## get group levels
  grpVar  <- obj$callArgs$groupVar
  grpLvls <- unique( obj$cuminc[[grpVar]] ) 

  ## check that 'refLvl' is one of the group levels found in the input object
  if ( !(refLvl %in% grpLvls) ) 
    stop("The value specified for 'refLvl' does not match any group levels found in input\n",
         "object '", as.list(match.call())$obj, "' \n\n")

  ## if cmpLvl is specified, check that the value(s) specified are among the group levels 
  if ( !is.null( cmpLvl ) ) {
    notFound <- !( cmpLvl %in% grpLvls ) 
    if ( any(notFound) ) 
      stop("The following levels specified via argument 'cmpLvl' were not found in the input:\n",
            paste( cmpLvl[ notFound ], collapse=","), ".\n\n" )
  }


  ## if cmpLvl wasn't specified, assign it's default value (all non-reference levels
  ## of the group variable )
  if ( is.null(cmpLvl) ) {
    cmpLvl <- grpLvls[ grpLvls != refLvl ]
  }

  # extract data.frame with cuminc estimates from input object
  cuminc <- obj$cuminc

  ## make vector of columns of cuminc to keep:

  ## get name of strata-variable, if one was used
  strataVar <- obj$callArgs$strataVar

  ## if a strata-var was used but is not in cuminc, set strataVar to NULL
  if ( !is.null(strataVar) ) {
    if ( !(strataVar %in% names(cuminc)) ) strataVar <- NULL
  }
 

  ## create var eventType as either 'eventType' or NULL
  eventType <- if ( inherits(obj, "crCumInc" ) ) "eventType" else NULL

  keepCols  <- c( strataVar, eventType, grpVar, "time", "cuminc", "se.cuminc" )
  matchCols <- c( strataVar, eventType, "time" )

  # then separate the reference level data from the rest  
  ref <- cuminc[ which( cuminc[[grpVar]] == refLvl ), keepCols ]
  cmp <- cuminc[ which( cuminc[[grpVar]] %in% cmpLvl ), keepCols ]

  ## Now we can merge the two data.frames
  joinDF <- merge(cmp, ref, by= matchCols,
                  all=TRUE, no.dups=TRUE, suffixes=c(".cmp",".ref") )

  alpha <- 1 - confLevel
  CIR <- within( joinDF, {

           cmpLvl = get( paste0( grpVar, ".cmp") )
           refLvl = get( paste0( grpVar, ".ref") )
           comparison = paste0( cmpLvl, " vs. ", refLvl)

           ## indicator of values where one cuminc is 0
           cumInc0 <- ( cuminc.cmp == 0  |  cuminc.ref == 0 )

           ## compute the quantities, ignoring presence of division-by-zero
           ##  and log-of-zero which will generate NaN values
           logCIR    <- log( cuminc.cmp ) - log( cuminc.ref )
           se.logCIR <- sqrt( (se.cuminc.cmp/cuminc.cmp)^2 + 
                              (se.cuminc.ref/cuminc.ref)^2  )

           ## replace NaN caused by zeros in above step
           logCIR[ cumInc0 ]    <- NA
           se.logCIR[ cumInc0 ] <- NA

           lo.logCIR <- logCIR - qnorm(1 - alpha/2)*se.logCIR
           up.logCIR <- logCIR + qnorm(1 - alpha/2)*se.logCIR


           ## if logCIR is defined, then exponentiate it to get CIR.  If logCIR is not 
           ## defined, CIR still may be - if it equals 0. We implement this below
           CIR <- ifelse( !is.na(logCIR), exp( logCIR ), ifelse( cuminc.ref > 0, 0L, NA ) )
           lo.CIR <- exp( lo.logCIR )
           up.CIR <- exp( up.logCIR )

           eff    <- 1 - CIR
           lo.eff <- 1 - up.CIR
           up.eff <- 1 - lo.CIR

           rm( cumInc0 )
       } )


  ## reorder components
  ordCols <- c(strataVar, eventType, "comparison", "cmpLvl", "refLvl", "time",
               "CIR", "lo.CIR", "up.CIR", "eff", "lo.eff", "up.eff", 
               "logCIR", "se.logCIR", "lo.logCIR", "up.logCIR", "cuminc.cmp", "cuminc.ref")
  CIR <- CIR[, ordCols]

  ## reorder rows
  row.ordr <- do.call("order", args=CIR[c(strataVar, "comparison", eventType, "time")] )
  CIR <- CIR[ row.ordr, ]


  ## slap CIR onto the object we were given, as a new component and then return it
  obj$CIR <- CIR

  ## Put results from the last timepoint (from each comparison) into a separate component,
  ## They are what will be shown as "the" estimates of VE (PE)
  effKeepVars <- c(strataVar, eventType, "comparison", "cmpLvl", "refLvl", 
                   "time", "eff", "lo.eff", "up.eff")

  # !! Will this keep the last time for each strata?
  obj$eff <- cbind( CIR[ CIR$time == max(CIR$time), effKeepVars ], confLevel = confLevel )


  # ---------------------------------
  # Derive Wald test of nullHypEff:
  # ---------------------------------
  # (1) get last row of CIR (per-strata) 
  # (2) extract components needed
  # (3) compute test stat (on log CIR scale, since that's what we have se computed for)
  # --------------------------------------------------------------------------------

  # get last row
  lastCIR   <- CIR[ CIR$time == max(CIR$time), ]

  ## null hypotheses on efficacy scale
  nullHypEff <- if ( is.null(nullHypEff)) 0L 
                else nullHypEff

  ## null hypotheses translated to the logCIR scale
  hypLogCIR <- log(1-nullHypEff) 

  # compute test statistics for two-sided tests of H0: Eff = nullHypEff
  # equivalent to H0: logCIR = log(1-nullHypEff) 
  testStat <- ( (hypLogCIR - lastCIR$logCIR)/lastCIR$se.logCIR )^2

  # compute the pvalues for the test-statistics
  pvalue <- 1 - pchisq( testStat, df=1 )

  tests <- data.frame(
              test = "Wald",
              nullHypEff = nullHypEff,
              nullHypHR  = 1 - nullHypEff,
              testStatistic = testStat,
              df = 1,
              pvalue = pvalue,
              stringsAsFactors = FALSE
           )

  ## add 'tests' onto output object
  obj$tests <- tests

  ## finally, store the call arguments
  obj$argsCIR <- list(refLvl = refLvl, cmpLvl = cmpLvl)

  ## add class 'CIR' onto the object
  class( obj ) <- c("CIR", class(obj))

  return( obj )
}





## Get cuminc curves for each event type, when there are competing risks
cmpRisks <- 
  function(data, futimeVar, eventVar, groupVar=NULL,  strataVar= NULL, 
           censor=list(minAtRisk= NULL, time= NULL), confLevel=0.95 )
{

##  'futimeVar' is a character string naming the variable in 'data' that contains
##      participant follow-up times
##
##  'eventVar' is a character string naming the variable in 'data' that contains
##      information on whether an event occurred.  If the named variable is not 
##      a factor, it will be coerced to one.  The FIRST level of the factor 
##      (i.e. levels( data[[eventVar]] )[1] ) will be assumed to represent 
##      "censoring" (i.e. "no event"), and the other levels will represent different
##      types of events.  It is highly recommended that *YOU* transform eventVar into
##      a factor before running cmpRisks() so you have control over the ordering of
##      the factor levels - else you may not get correct results.
##
##  'groupVar' is a character string giving the name of a variable in 'data'. 
##     The variable must be categorical and should identify some characteristic of interest.
##     Typically this would be something like treatment group, but could be anything.
##
##  'strataVar' is a character string giving the name of the variable in 'data' to be
##     stratified on.  If you have multiple variables you want to stratify on (e.g. sex 
##     and dosage) you can do this by first creating a variable in your dataset that 
##     jointly specifies the levels of all those variables and pass that variable name
##     via 'strataVar'.  You can create such a variable using e.g. interaction(). 
##     Note that, 'stratification' in this function will simply allow separate estimation  
##     of the cumulative incidence functions for different levels of these variables,
##     there is currently no method for combining the estimates from different stratas
##     to create a "stratified-estimate".
##
##  'confLevel' gives the desired confidence level for the confidence intervals for: 
##     cumhaz, cuminc, and cir/PE if requested 
##
##  'censor' allows a way to specify if/how censoring should be done.  If 'censor' is 
##     specified, it must be a list with components named either 'minAtRisk' or 'time'.
##     If 'minAtRisk' = m is specified, then that instructs that censoring be done to 
##     ensure that the number of participants at risk at all timepoints is >= m.  This is
##     implemented by finding the latest timepoint such that all group x strata subcohorts
##     satisfy the 'minimum-at-risk' criteria at that timepoint - that timepoint will be 
##     returned as 'censorTime'.  If 'time' = t is specified, then the data will be censored
##     at time 't'.  If *both* 'minAtRisk' and 'time'=t are specified, the time that satisfies 
##     the 'minAtRisk' criteria will be determined, and then the *smaller* of that time and
##     time 't' will be used for censoring (and *that* value returned via 'censorTime')


  ##  'eps' (epsilon) = a small number used to avoid issues of floating point arithmetic
  eps <- sqrt( .Machine$double.eps )


  ## check all "Var" inputs - must be character vectors and be in 'data'
  if ( !is.null(groupVar) && ( length(groupVar) > 1 || !is.character(groupVar) ) )
    stop("Error: argument 'groupVar' must be a character vector of length 1\n\n")

  if ( !is.null(strataVar) && ( length(strataVar) > 1 || !is.character(strataVar) ) )
    stop("Error: argument 'strataVar' must be a character vector of length 1\n\n")

  if ( length(futimeVar) > 1 || !is.character(futimeVar) )
    stop("Error: argument 'futimeVar' must be a character vector of length 1\n\n")

  if ( length(eventVar) > 1 || !is.character(eventVar) )
    stop("Error: argument 'eventVar' must be a character vector of length 1\n\n")

  if (!is.null(groupVar) && !groupVar %in% names(data) )
    stop("The variable ", groupVar, " was not found in dataset ", data, "\n\n")

  if (! futimeVar %in% names(data) )
    stop("The variable ", futimeVar, " was not found in dataset ", data, "\n\n")

  if (! eventVar %in% names(data) )
    stop("The variable ", eventVar, " was not found in dataset ", data, "\n\n")

  if (!is.null(strataVar) && !strataVar %in% names(data) )
    stop("The variable ", strataVar, " was not found in dataset ", data, "\n\n")


  if ( !missing(censor) ) {
    if ( !is.list(censor) )
      stop("Error: argument 'censor' must be a list\n\n")

    censor.components <- c('minAtRisk', 'time')
    if ( !all( names(censor)  %in% censor.components ) )
      stop("Error: 'censor' may only have compents with names among: ",
           paste(censor.components, collapse=", "),"\n\n")

    if (!is.null(censor$minAtRisk)) {
       ## ensure minAtRisk is an integer. We allow a non-integer to be specified,
       ## but if one is, it's converted here to the next higher integer
       censor$minAtRisk <- ceiling( censor$minAtRisk )
    }
  } else {
    ## if censor wasn't specified, set it to NULL
    censor <- NULL
  }


  ## Now that we're done checking inputs, let's define our output object structure and 
  ## fill in some pieces
  obj <- list(
              callArgs = list(
                  futimeVar = futimeVar,
                  eventVar  = eventVar,
                  groupVar  = groupVar,
                  strataVar = strataVar,
                  censor = censor,
                  confLevel = confLevel
              )[ c(TRUE, TRUE, !is.null(groupVar), !is.null(strataVar), !missing(censor), TRUE ) ],


              ## indicator of whether any censoring was *PERFORMED* on the data
              ## (as opposed to having simply specifyied censoring arguments - since
              ##  it's possible that they'd have no effect)
              censored = FALSE,

              ## summary info on censoring - will be filled in if censoring was ATTEMPTED
              censorInfo = list(
                  censorTime = NA,
                  nRecsCensored = NA,
                  nEventsCensored = NA
                ), 

              ## the data used in estimation is returned here, will be censored if 
              ## output component 'censored' is TRUE
              data = NULL,

              ## data.frame of cumulative incidence and cumulative hazard estimates,
              ## standard errors and confidence intervals, possibly for separate strata, 
              ## e.g. if 'strataVar' was specified but a stratified estimator was not
              ## requested 
              cuminc = NULL 
            )

  ## put a class on the object (the "na" = "nelson-aalen")
  class(obj) <- c("crCumInc", class(obj)) 


  ## post-process 'censor'
  ## ---------------------
  if ( !is.null(censor) ) {
    ## if it was specified, remove NULL components (which could be both of them!)
    null.components <- sapply(censor, is.null)

    ## if they're all NULL, set censor to NULL
    if ( all( null.components ) ) {
       censor <- NULL
    } else {
      ## Keep  non-NULL components
      censor <- censor[ !null.components ]
    }
  }

  ## subset the data to the variables we need, prior to modifying
  d <- data[, c(futimeVar, eventVar, groupVar, strataVar) ]


  ## is 'eventVar' a factor?  If not, convert to a factor (since survival will
  ## anyways) and issue a NOTE that this was done and what the censoring level will be
  if ( !is.factor( d[[ eventVar ]] ) ) {
    cat("\n\n", "WARNING:\n", 
       "The competing risks analysis requires that the event-variable be a factor, \n",
       "the one you've provided is not and so must be coerced to one.  This is potentially\n",
       "problematic since the first factor levels is assumed represent 'no event'/censoring \n",
       "and the ordering of the levels is affected by the user's environment.  So different \n",
       "may get different results, or the same user may get different results on different \n",
       "machines.  Bottom line: you should convert this variable to a factor BEFORE passing \n",
       "the data to this function - and explicitly specify the level ordering when creating",
       "it \n\n\n", sep="")

    cat("Converting ", eventVar, " to a factor...\n\n", sep="")

    d[[ eventVar ]] <- as.factor( d[[ eventVar ]] )

    ## get first level of the factor
    first.level <- levels( d[[ eventVar ]])[1]

    cat("NOTE: the first level of the factor ", eventVar, " is: ", first.level, "- it \n",
        "is assumed to represent censored outcomes.  If this is not correct, please \n",
        "convert ", eventVar, " to a factor before runnning this function, and explicitly \n",
        "specify the level ordering when doing so", "\n\n\n", sep="")
  }

  if (is.null(groupVar) && is.null(strataVar) ) {

    ## 'Strata' is created for internal use
    Strata <- FALSE
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ 1")

  } else {
    ## create our own temporary interaction variable, combining 'groupVar' and
    ## 'strataVar'.  This will allow us to unwind some of the changes that will
    ## be made to the levels in the output of survfit.
    ## The code below works whether or not 'strataVar' is specified
    tmpStrataVar <- "my_tmp_Strata"
    d[[tmpStrataVar]] <- interaction(d[,c(groupVar,strataVar)], lex.order=TRUE)

    ## store the mapping of the new strata var. to its component ones
    strataMapping <- unique(d[, c(tmpStrataVar, groupVar, strataVar)] )

    Strata <- TRUE
    formula.char <- paste0( "Surv(", futimeVar, ",", eventVar, ") ~ ", tmpStrataVar)
  }


  ## Now we're ready to do censoring, if any was requested
  if ( !is.null(censor) ) { 
    if ( is.null( censor$minAtRisk ) ) {
      ## if minAtRisk not given, then we have only censor$time to work with - EZ stuff!
      censorTime <- censor$time
    } else {     
      ## determine the censorTime corresponding to minAtRisk 
      minAR <- censor$minAtRisk

      ## use 'minAR' to determine the time to censor at
      if (!Strata) {
        # if no groups/strata, use all data 
        censorTime <- sort( d[[ futimeVar ]], decreasing=TRUE )[ minAR ] 

      } else {

        # else do by strata
        censTimes.strata <- 
            by(data = d, INDICES = d[[tmpStrataVar]], 
               FUN = function(x,n) sort(x[[ futimeVar ]], decreasing=TRUE )[n],
               n = minAR )

        censorTime <- min( as.vector(censTimes.strata), na.rm=TRUE )
      }

      ## if censor$time was *also* specified, then the minimum of it and the time just computed
      if ( !is.null( censor$time ) ) 
        censorTime <- min( censorTime, censor$time, na.rm=TRUE )
    }

    ## Now we actually *do* the censoring, censoring the data at 'censorTime'
    d <- censorData( 
            data = d,
            censorTime= censorTime, 
            futimeVar= futimeVar, 
            eventVar= eventVar,
            censInd = TRUE ) 

    ## set censoring indicator in output object
    obj$censored <- ifelse( sum(d$censored) > 0, TRUE, FALSE ) 

    ## set the censoring time that was used (regardless of whether it had an effect)
    obj$censorInfo$censorTime <- censorTime

    ## if censoring was done, modify data, compute some info and store
    if ( obj$censored ) {

      ## also attach the original data futime and event indicator columns to the
      ## censored dataset, storing them in variables named by adding a "orig_" prefix
      ## to the original variable name
      d[[ paste0("orig_", futimeVar) ]] <- data[[ futimeVar ]]
      d[[ paste0("orig_", eventVar) ]]  <- data[[ eventVar ]]

      censLvl <- levels( d[[eventVar]] )[1]
      obj$censorInfo$nRecsCensored   <- sum(d$censored)
      obj$censorInfo$nEventsCensored <- 
          sum( d[[ paste0("orig_", eventVar) ]] != d[[ eventVar ]]  )

    } else {

      ## no censoring done
      obj$censorInfo$nRecsCensored   <- 0
      obj$censorInfo$nEventsCensored <- 0

      ## remove censoring indicator from data
      d$censored <- NULL
    } 
  }
  ## store data into output object (excluding the temporary strata var we created, if present)
  obj$data <- d[ , names(d) != ifelse(Strata, tmpStrataVar,"") ] 



  ## load survival package
  library(survival)

  ## coerce our character string 'formula.char' to an actual "formula"
  SurvFormula <- as.formula( formula.char )


  ## identify all "event times", which are the levels of eventVar that not equal to
  ## the factor's first level (first level represents 'censoring'/'no event')
  censLvl <- levels( d[[eventVar]] )[1]
  eventInd <- d[[ eventVar ]] != censLvl
  all.uniq.event.times <- c(0, sort( unique( d[eventInd, futimeVar]) ) )

  summ <- summary( survfit( SurvFormula, data=d, conf.int=confLevel),
                   extend = TRUE, times = all.uniq.event.times )

  ## number of types of events
  nEvents <- nlevels( d[[eventVar]] ) - 1

  dfList <- vector("list", length=nEvents)
  for ( i in 1:nEvents ) {

    summ.i <- data.frame(
              eventType = summ$states[i],
              time      = summ$time,
              n.risk    = summ$n.risk[, nEvents+1],
              n.event   = summ$n.event[, i],
              cuminc    = summ$pstate[, i],
              lo.cuminc = summ$lower[, i],
              up.cuminc = summ$upper[, i],
              se.cuminc = summ$std.err[, i],
              stringsAsFactors = FALSE )
    if ( Strata )
      summ.i$strata <- summ$strata

    dfList[[ i ]] <- summ.i
  }

  cumInc <- do.call(rbind, dfList )

  if ( Strata ) {
    ## fix the levels of 'strata':
    ## survfit alters the level names of 'strata', our code (below) changes them back
    levels( cumInc$strata ) <- sub( paste0(tmpStrataVar,"="), "", levels(cumInc$strata) )

    ## rename 'strata' before merging to add original variables back on (in case one
    ## of the original variables has that same name)
    names(cumInc)[ names(cumInc) == "strata" ] <- tmpStrataVar

    ## deconvolute the 'Strata' variable to extract 'group, and then drop the temp variable
    cumInc <- merge( strataMapping, cumInc ) #, by= get("tmpStrataVar") )

    ## re-sort the data, by tmpStrataVar and time
    cumInc <- cumInc[ order(cumInc[[tmpStrataVar]], cumInc$eventType, cumInc$time), ] 

    ## drop the tmp strata variable
    cumInc[[ tmpStrataVar ]] <- NULL
  } else {

    ## re-sort the data, by groupVar and time
    cumInc <- cumInc[ order(cumInc[[groupVar]], cumInc$eventType, cumInc$time), ] 
  }

  ## re-org the columns (and exclude ones we don't want to keep)

  ## Add 'eventType' variable here too
  ordCols <- c(groupVar,  strataVar, "eventType", "time", "n.risk", "n.event",
               "cuminc", "se.cuminc", "lo.cuminc", "up.cuminc" )
  cumInc <- cumInc[ , ordCols ]

  ## if a stratified estimator was NOT requested, then we're done here
  obj$cuminc <- cumInc

  return( obj )
}
