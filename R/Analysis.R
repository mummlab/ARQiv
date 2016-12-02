#############################################
# SIGNAL EXTRACTION
#############################################
#' @title Signal Extraction
#' @description This function reads .xml input files, extracts and isolates the signal from individual wells of a 96-well plate. This function is applied to both control and compound plates.
#'
#' @param xmlfile A .xml file provided by the plate reader containing raw data from regional scanns of each well.
#' @return Numerical 'signal' values (arbitrary units) for all individual wells which are applied to downstream analyses.
ParseXML <- function(xmlfile=NULL)
{
  xmlline <- readLines(xmlfile);
  # RETRIEVE THE XML FILE NAME
  exprm.name <- xmlline[grep("<ID>", xmlline)];
  exprm.name <- gsub(" ", "", exprm.name);
  exprm.name <- gsub("<ID>", "", exprm.name);
  exprm.name <- gsub("</ID>", "", exprm.name);
  # RETRIEVE THE SIGNAL OF EACH WELL
  all.well.signal <- NULL;
  well.pos <- grep("Well Pos", xmlline);
  for(i in 1:length(well.pos))
  {
    well.signal <- array(0, dim=11);
    # 0;2 SIGNAL
    well.signal[1] <- xmlline[well.pos[i]+1];
    well.signal[1] <- gsub(" ", "", well.signal[1]);
    well.signal[1] <- gsub("<MultipleMRW_Position=\"0;2\"LabelId=\"1\">", "", well.signal[1]);
    well.signal[1] <- as.numeric(gsub("</Multiple>", "", well.signal[1]));
    # 1;2 SIGNAL
    well.signal[2] <- xmlline[well.pos[i]+2];
    well.signal[2] <- gsub(" ", "", well.signal[2]);
    well.signal[2] <- gsub("<MultipleMRW_Position=\"1;2\"LabelId=\"1\">", "", well.signal[2]);
    well.signal[2] <- as.numeric(gsub("</Multiple>", "", well.signal[2]));
    # 2;2 SIGNAL
    well.signal[3] <- xmlline[well.pos[i]+3];
    well.signal[3] <- gsub(" ", "", well.signal[3]);
    well.signal[3] <- gsub("<MultipleMRW_Position=\"2;2\"LabelId=\"1\">", "", well.signal[3]);
    well.signal[3] <- as.numeric(gsub("</Multiple>", "", well.signal[3]));
    # 2;1 SIGNAL
    well.signal[4] <- xmlline[well.pos[i]+4];
    well.signal[4] <- gsub(" ", "", well.signal[4]);
    well.signal[4] <- gsub("<MultipleMRW_Position=\"2;1\"LabelId=\"1\">", "", well.signal[4]);
    well.signal[4] <- as.numeric(gsub("</Multiple>", "", well.signal[4]));
    # 1;1 SIGNAL
    well.signal[5] <- xmlline[well.pos[i]+5];
    well.signal[5] <- gsub(" ", "", well.signal[5]);
    well.signal[5] <- gsub("<MultipleMRW_Position=\"1;1\"LabelId=\"1\">", "", well.signal[5]);
    well.signal[5] <- as.numeric(gsub("</Multiple>", "", well.signal[5]));
    # 0;1 SIGNAL
    well.signal[6] <- xmlline[well.pos[i]+6];
    well.signal[6] <- gsub(" ", "", well.signal[6]);
    well.signal[6] <- gsub("<MultipleMRW_Position=\"0;1\"LabelId=\"1\">", "", well.signal[6]);
    well.signal[6] <- as.numeric(gsub("</Multiple>", "", well.signal[6]));
    # 0;0 SIGNAL
    well.signal[7] <- xmlline[well.pos[i]+7];
    well.signal[7] <- gsub(" ", "", well.signal[7]);
    well.signal[7] <- gsub("<MultipleMRW_Position=\"0;0\"LabelId=\"1\">", "", well.signal[7]);
    well.signal[7] <- as.numeric(gsub("</Multiple>", "", well.signal[7]));
    # 1;0 SIGNAL
    well.signal[8] <- xmlline[well.pos[i]+8];
    well.signal[8] <- gsub(" ", "", well.signal[8]);
    well.signal[8] <- gsub("<MultipleMRW_Position=\"1;0\"LabelId=\"1\">", "", well.signal[8]);
    well.signal[8] <- as.numeric(gsub("</Multiple>", "", well.signal[8]));
    # 2;0 SIGNAL
    well.signal[9] <- xmlline[well.pos[i]+9];
    well.signal[9] <- gsub(" ", "", well.signal[9]);
    well.signal[9] <- gsub("<MultipleMRW_Position=\"2;0\"LabelId=\"1\">", "", well.signal[9]);
    well.signal[9] <- as.numeric(gsub("</Multiple>", "", well.signal[9]));
    # Mean SIGNAL
    well.signal[10] <- xmlline[well.pos[i]+10];
    well.signal[10] <- gsub(" ", "", well.signal[10]);
    well.signal[10] <- gsub("<MultipleMRW_Position=\"Mean\"LabelId=\"1\">", "", well.signal[10]);
    well.signal[10] <- as.numeric(gsub("</Multiple>", "", well.signal[10]));
    # StDev SIGNAL
    well.signal[11] <- xmlline[well.pos[i]+11];
    well.signal[11] <- gsub(" ", "", well.signal[11]);
    well.signal[11] <- gsub("<MultipleMRW_Position=\"StDev\"LabelId=\"1\">", "", well.signal[11]);
    well.signal[11] <- as.numeric(gsub("</Multiple>", "", well.signal[11]));
    well.signal <- as.numeric(well.signal);
    all.well.signal <- rbind(all.well.signal, well.signal);
  }
  rownames(all.well.signal) <- c(paste0("A", 1:12), paste0("B", 12:1), paste0("C", 1:12), paste0("D", 12:1), paste0("E", 1:12), paste0("F", 12:1), paste0("G", 1:12), paste0("H", 12:1));
  all.well.signal <- all.well.signal[c(1:12, 24:13, 25:36, 48:37, 49:60, 72:61, 73:84, 96:85), ];
  colnames(all.well.signal) <- c("S02", "S12", "S22", "S21", "S11", "S01", "S00", "S10", "S20", "Mean", "StDev");
  exprm <- list();
  exprm.name <- gsub("/", "-", exprm.name)
  exprm$name <- exprm.name;
  exprm$signal <- all.well.signal;
  return(exprm);
}

#############################################
# BACKGROUND SIGNAL CALCULATION
#############################################
#' @title Background Signal Calculation
#' @description This function calculates background from non-transgenic larvae. They identical to their transgenic counterparts with regard to pigmentation (e.g., Roy Orbison pigmentation mutants).
#'
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @return bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @export


SIGNAL <- function(Roy.file=NULL){
  if(is.null(Roy.file)) stop("Please input your Roy.file")
  Roy <- ParseXML(xmlfile=Roy.file);
  Roy.signal <- Roy$signal[,1:9];
  Roy.mean <- mean(apply(Roy.signal, 1, max));
  Roy.sd <- sd(apply(Roy.signal, 1, max));

  #output is minimal signal value
  bg.ratio <- Roy.mean+3*Roy.sd;
  return(bg.ratio)
}




#############################################
# SAMPLE SIZE DETERMINATION
#############################################
#' @title Sample Size Determination
#' @description This function calculates the minimum sample size required for
#' your assay with varying paremeters, such as simulated effect size
#' and TypeI/TypeII error rates. This can be completed using .XML or .CSV file formats.
#'
#' @param Pos.file A .xml file derived from your positive control
#'  condition, i.e., that with the highest value.
#' @param Neg.file A .xml file derived from your negative
#' control condition, i.e., that with the lowest value.
#' @param Excel.file A .xlsx file derived manually containing
#' single columns of 'signal' from both positive and negative controls.
#' @param Roy.file A .xml file of raw signal values from a
#' 96-well plate of non-transgenic Roy larvae.
#' @param effect.size A simulated effect size for virtual
#' experiments (e.g., 0.75 is a compound that simulates 75 percents efficacy when compared to the positive control condition).
#' @param bg.ratio Numerical minimal 'signal' cutoff value
#' (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' Default value is 3946.
#' @param dic.output Output directory.Default is current working directory.
#' @return SampleSize.csv A .csv file of minimum sample sizes required
#' to observe simulated effect sizes based on varied TypeI/TypeII errors.
#' @import utils
#' @export
SAMPLE <- function(Pos.file=NULL, Neg.file=NULL,   Roy.file=NULL,Excel.file=NULL, effect.size=NULL,bg.ratio=NULL,dic.output=NULL)
{
  if(!is.null(Roy.file)){
    bg.ratio <- SIGNAL(Roy.file)
  }else{
    if(is.null(bg.ratio)) bg.ratio <- 3946
  }
  if(is.null(Excel.file)){
    if(is.null(Pos.file)) Pos.file <- system.file("extdata","Positive.Control.xml",package = "ARQiv")
    if(is.null(Neg.file)) Neg.file <- system.file("extdata","Negative.Control.xml",package = "ARQiv")


    # POSITIVE CONTROL SIGNAL EXTRACTION
    Pos.XML <- ParseXML(xmlfile=Pos.file);
    Pos.signal <- as.matrix(Pos.XML$signal[ ,1:9]);
    Pmax.signal <- apply(Pos.signal, 1, max);
    Pos.signal[Pos.signal<=bg.ratio] <- 0;
    Pwell.signal <- apply(Pos.signal, 1, sum);
    Pwell.signal <- as.numeric(apply(cbind(Pwell.signal, Pmax.signal), 1, max));
    s.pos <- Pwell.signal;

    # NEGATIVE CONTROL SIGNAL EXTRACTION
    Neg.XML <- ParseXML(xmlfile=Neg.file);
    Neg.signal <- as.matrix(Neg.XML$signal[ ,1:9]);
    Nmax.signal <- apply(Neg.signal, 1, max);
    Neg.signal[Neg.signal<=bg.ratio] <- 0;
    Nwell.signal <- apply(Neg.signal, 1, sum);
    Nwell.signal <- as.numeric(apply(cbind(Nwell.signal, Nmax.signal), 1, max));
    s.neg <- Nwell.signal;
  }else{
    excel.signal <- read.csv(Excel.file);
    s.pos <- excel.signal$Pos[!is.na(excel.signal$Pos)];
    s.neg <- excel.signal$Neg[!is.na(excel.signal$Neg)];
  }

  # SAMPLE SIZE DETEMINATION FOR PERCENT EFFECT
  if(is.null(effect.size)){
    Decrease.ratio <- c("100%", "75%", "50%", "25%")
  }else{
    effect.size <- as.vector(effect.size)
    effect.size <- pmin(effect.size,1)
    effect.size <- pmax(effect.size,0)
    Decrease.ratio <- paste0(effect.size*100,"%")
  }

  TypeI.error <- c(0.05, 0.01, 0.001);
  TypeII.error <- c(0.2, 0.05, 0.01, 0.001);
  Z.a <- c(1.96, 2.58, 3.29);
  Z.b <- c(0.84, 1.64, 2.33, 3.09);
  len.Za <- length(Z.a);
  len.Zb <- length(Z.b);
  if(is.null(effect.size)){
    fold75 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.75+mean(s.neg));
    fold50 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.5+mean(s.neg));
    fold25 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.25+mean(s.neg));
    Decrease.fold <- c(1, fold75, fold50, fold25)
  }else{
    Decrease.fold <- (mean(s.pos)/((mean(s.pos)-mean(s.neg))*effect.size+mean(s.neg)))
  }
  TypeI.error <- rep(TypeI.error, rep(c(len.Zb*length(Decrease.fold)), len.Za));
  TypeII.error <- rep(rep(TypeII.error, rep(length(Decrease.fold),len.Zb)), len.Za);
  Z.a <- rep(Z.a, rep(c(len.Zb*length(Decrease.fold)), len.Za));
  Z.b <- rep(rep(Z.b,rep(length(Decrease.fold),len.Zb)), len.Za);

  smpl.size <- array(0, dim=length(Decrease.ratio));
  Decrease.ratio <- rep(Decrease.ratio, time=len.Za*len.Zb);
  Decrease.fold <- rep(Decrease.fold, time=len.Za*len.Zb);
  for(i in 1:length(Decrease.ratio))
  {
    P.decrs <- s.pos/Decrease.fold[i];
    smpl.size[i] <- 2*max(sd(P.decrs), sd(s.neg))^2*(Z.a[i]+Z.b[i])^2/(mean(P.decrs)-mean(s.neg))^2*1.15
    smpl.size[i] <- round(smpl.size[i],2)
  }
  log.smpl.size <- array(0, dim=length(Decrease.ratio));
  for(i in 1:length(Decrease.ratio))
  {
    P.decrs <- s.pos/Decrease.fold[i];
    log.smpl.size[i] <- 2*max(sd(log(P.decrs,base=2)), sd(log(s.neg, base=2)))^2*(Z.a[i]+Z.b[i])^2/(mean(log(P.decrs, base=2))-mean(log(s.neg, base=2)))^2*1.15
    log.smpl.size[i] <- round(log.smpl.size[i],2)
  }

  sample.size <- data.frame(Decrease.ratio, Decrease.fold, TypeI.error, TypeII.error, Z.a, Z.b, smpl.size, log.smpl.size);
  colnames(sample.size) <- c("Decrease Ratio", "Decrease Fold", "TypeI error", "TypeII error", "Z(Alpha)", "Z(Beta)", "Sample Size", "Sample Size(log)");
  if(is.null(dic.output)) dic.output <- getwd()
  write.csv(sample.size, file=paste0(dic.output,"/SampleSize.csv"), row.names=FALSE);
}


#############################################
# ZFACTOR
#############################################
#' @title Zfactor Calculation
#' @description This function calculates SSMD_EST of your chosen assay
#' parameters, using both positive and negative control data. This can be
#'  completed using .XML or .CSV file formats.
#'
#' @param Pos.file A .xml file derived from your positive control condition, i.e., that with the highest value.
#' @param Neg.file A .xml file derived from your negative control condition, i.e., that with the lowest value.
#' @param Excel.file A .xlsx file derived manually containing single columns of 'signal' from both positive and negative controls.
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @param effect.size A simulated effect size for virtual experiments (e.g., 0.75 is a compound that simulates 75 percents efficacy when compared to the positive control condition).
#' @param bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @param sampling.num The number of invidiual wells you will sample from a 96-well plate.
#' @param times The number of iterative cycles of the chosen sampling number.
#' @param dic.output Output directory.Default is current working directory.
#' @return a .csv file with SSMD_EST estimated via random sampling (i.e., a compound with 75% efficacy when compared to the positive control will have an SSMD score of X).
#' @export

ZFACTOR <- function(Pos.file =NULL, Roy.file=NULL, bg.ratio=NULL,Neg.file=NULL,Excel.file=NULL,dic.output=NULL)
{

  if(!is.null(Roy.file)){
    bg.ratio <- SIGNAL(Roy.file)
  }else{
    if(is.null(bg.ratio)) bg.ratio <- 3946
  }

  if(is.null(Pos.file)|is.null(Neg.file)){
    if(is.null(Excel.file)) Excel.file <- system.file("extdata","Excel.Pos.Neg.Control.csv",package = "ARQiv")
    excel.signal <- read.csv(Excel.file);
    s.pos <- excel.signal$Pos[!is.na(excel.signal$Pos)];
    s.neg <- excel.signal$Neg[!is.na(excel.signal$Neg)];

  }else{
    if(is.null(Pos.file)) Pos.file <- system.file("extdata","Positive.Control.xml",package = "ARQiv")
    if(is.null(Neg.file)) Neg.file <- system.file("extdata","Negative.Control.xml",package = "ARQiv")

    # POSITIVE CONTROL SIGNAL EXTRACTION
    Pos.XML <- ParseXML(xmlfile=Pos.file);
    Pos.signal <- as.matrix(Pos.XML$signal[ ,1:9]);
    Pmax.signal <- apply(Pos.signal, 1, max);
    Pos.signal[Pos.signal<=bg.ratio] <- 0;
    Pwell.signal <- apply(Pos.signal, 1, sum);
    Pwell.signal <- as.numeric(apply(cbind(Pwell.signal, Pmax.signal), 1, max));
    s.pos <- Pwell.signal;

    # NEGATIVE CONTROL SIGNAL EXTRACTION
    Neg.XML <- ParseXML(xmlfile=Neg.file);
    Neg.signal <- as.matrix(Neg.XML$signal[ ,1:9]);
    Nmax.signal <- apply(Neg.signal, 1, max);
    Neg.signal[Neg.signal<=bg.ratio] <- 0;
    Nwell.signal <- apply(Neg.signal, 1, sum);
    Nwell.signal <- as.numeric(apply(cbind(Nwell.signal, Nmax.signal), 1, max));
    s.neg <- Nwell.signal;
  }


  P.signal <- s.pos
  N.signal <- s.neg

  P.signal.o <- P.signal
  z.factor <- 1-3*(sd(P.signal)+sd(N.signal))/abs(mean(P.signal)-mean(N.signal))

  return(round(z.factor,2))
}


#############################################
# SSMD_QC ESTIMATION
#############################################
#' @title SSMD_QC Calculation
#' @description This function calculates SSMD_EST of your chosen
#' assay parameters, using both positive and negative control
#'  data. This can be completed using .XML or .CSV file formats.
#'
#' @param Pos.file A .xml file derived from your positive control condition, i.e., that with the highest value.
#' @param Neg.file A .xml file derived from your negative control condition, i.e., that with the lowest value.
#' @param Excel.file A .xlsx file derived manually containing single columns of 'signal' from both positive and negative controls.
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @param effect.size A simulated effect size for virtual experiments (e.g., 0.75 is a compound that simulates 75 percents efficacy when compared to the positive control condition).
#' @param bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @param sampling.num The number of invidiual wells you will sample from a 96-well plate.
#' @param times The number of iterative cycles of the chosen sampling number.
#' @param dic.output Output directory.Default is current working directory.
#' @return A .csv file with SSMD_EST estimated via random sampling (i.e., a compound with 75% efficacy when compared to the positive control will have an SSMD score of X).
#' @export

SSMD.QC <- function(Pos.file =NULL, Roy.file=NULL, bg.ratio=NULL,Neg.file=NULL,Excel.file=NULL,dic.output=NULL)
{

  if(!is.null(Roy.file)){
    bg.ratio <- SIGNAL(Roy.file)
  }else{
    if(is.null(bg.ratio)) bg.ratio <- 3946
  }

  if(is.null(Pos.file)|is.null(Neg.file)){

    if(is.null(Excel.file)) Excel.file <- system.file("extdata","Excel.Pos.Neg.Control.csv",package = "ARQiv")
    excel.signal <- read.csv(Excel.file);
    s.pos <- excel.signal$Pos[!is.na(excel.signal$Pos)];
    s.neg <- excel.signal$Neg[!is.na(excel.signal$Neg)];

  }else{
    if(is.null(Pos.file)) Pos.file <- system.file("extdata","Positive.Control.xml",package = "ARQiv")
    if(is.null(Neg.file)) Neg.file <- system.file("extdata","Negative.Control.xml",package = "ARQiv")


    # POSITIVE CONTROL SIGNAL EXTRACTION
    Pos.XML <- ParseXML(xmlfile=Pos.file);
    Pos.signal <- as.matrix(Pos.XML$signal[ ,1:9]);
    Pmax.signal <- apply(Pos.signal, 1, max);
    Pos.signal[Pos.signal<=bg.ratio] <- 0;
    Pwell.signal <- apply(Pos.signal, 1, sum);
    Pwell.signal <- as.numeric(apply(cbind(Pwell.signal, Pmax.signal), 1, max));
    s.pos <- Pwell.signal;

    # NEGATIVE CONTROL SIGNAL EXTRACTION
    Neg.XML <- ParseXML(xmlfile=Neg.file);
    Neg.signal <- as.matrix(Neg.XML$signal[ ,1:9]);
    Nmax.signal <- apply(Neg.signal, 1, max);
    Neg.signal[Neg.signal<=bg.ratio] <- 0;
    Nwell.signal <- apply(Neg.signal, 1, sum);
    Nwell.signal <- as.numeric(apply(cbind(Nwell.signal, Nmax.signal), 1, max));
    s.neg <- Nwell.signal;
  }

  P.signal <- s.pos
  N.signal <- s.neg


  ssmd.qc <- (median(P.signal)-median(N.signal))/(1.4826)/sqrt(mad(P.signal,constant = 1)^2+mad(N.signal,constant = 1)^2)
  return(round(ssmd.qc,2))
}

#############################################
# SSMD_EST ESTIMATION
#############################################
#' @title SSMD Estimation With Random Sampling
#' @description This function calculates SSMD_EST of your chosen assay parameters,
#' using both positive and negative control data. This can be completed using .XML or .CSV file formats.
#'
#' @param Pos.file A .xml file derived from your positive control condition, i.e., that with the highest value.
#' @param Neg.file A .xml file derived from your negative control condition, i.e., that with the lowest value.
#' @param Excel.file A .xlsx file derived manually containing single columns of 'signal' from both positive and negative controls.
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @param effect.size A simulated effect size for virtual experiments (e.g., 0.75 is a compound that simulates 75 percents efficacy when compared to the positive control condition).
#' @param bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @param sampling.num The number of invidiual wells you will sample from a 96-well plate.
#' @param times The number of iterative cycles of the chosen sampling number.
#' @param dic.output Output directory.Default is current working directory.
#' @return a .csv file with SSMD_EST estimated via random sampling (i.e., a compound with 75 percents efficacy when compared to the positive control will have an SSMD score of X).
#' @export

SSMD.EST <- function(Pos.file =NULL, Roy.file=NULL, bg.ratio=NULL,Neg.file=NULL, Excel.file=NULL, effect.size=NULL, sampling.num=16,times=10000,dic.output=NULL)
{
  if(!is.null(Roy.file)){
    bg.ratio <- SIGNAL(Roy.file)
  }else{
    if(is.null(bg.ratio)) bg.ratio <- 3946
  }
  if(is.null(Pos.file)|is.null(Neg.file)){

    if(is.null(Excel.file)) Excel.file <- system.file("extdata","Excel.Pos.Neg.Control.csv",package = "ARQiv")
    excel.signal <- read.csv(Excel.file);
    s.pos <- excel.signal$Pos[!is.na(excel.signal$Pos)];
    s.neg <- excel.signal$Neg[!is.na(excel.signal$Neg)];

  }else{
    if(is.null(Pos.file)) Pos.file <- system.file("extdata","Positive.Control.xml",package = "ARQiv")
    if(is.null(Neg.file)) Neg.file <- system.file("extdata","Negative.Control.xml",package = "ARQiv")

    # POSITIVE CONTROL SIGNAL EXTRACTION
    Pos.XML <- ParseXML(xmlfile=Pos.file);
    Pos.signal <- as.matrix(Pos.XML$signal[ ,1:9]);
    Pmax.signal <- apply(Pos.signal, 1, max);
    Pos.signal[Pos.signal<=bg.ratio] <- 0;
    Pwell.signal <- apply(Pos.signal, 1, sum);
    Pwell.signal <- as.numeric(apply(cbind(Pwell.signal, Pmax.signal), 1, max));
    s.pos <- Pwell.signal;

    # NEGATIVE CONTROL SIGNAL EXTRACTION
    Neg.XML <- ParseXML(xmlfile=Neg.file);
    Neg.signal <- as.matrix(Neg.XML$signal[ ,1:9]);
    Nmax.signal <- apply(Neg.signal, 1, max);
    Neg.signal[Neg.signal<=bg.ratio] <- 0;
    Nwell.signal <- apply(Neg.signal, 1, sum);
    Nwell.signal <- as.numeric(apply(cbind(Nwell.signal, Nmax.signal), 1, max));
    s.neg <- Nwell.signal;
  }
  P.signal <- s.pos
  N.signal <- s.neg

  if(is.null(effect.size)){
    Decrease.ratio <- c("100%", "75%", "50%", "25%")
  }else{
    effect.size <- as.vector(effect.size)
    effect.size <- pmin(effect.size,1)
    effect.size <- pmax(effect.size,0)
    Decrease.ratio <- paste0(effect.size*100,"%")
  }

  if(is.null(effect.size)){
    fold75 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.75+mean(s.neg));
    fold50 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.5+mean(s.neg));
    fold25 <- mean(s.pos)/((mean(s.pos)-mean(s.neg))*0.25+mean(s.neg));
    Decrease.fold <- c(1, fold75, fold50, fold25)
  }else{
    Decrease.fold <- (mean(s.pos)/((mean(s.pos)-mean(s.neg))*effect.size+mean(s.neg)))
  }

  P.signal.o <- P.signal
  sampling.num <- round(as.vector(sampling.num))[1]
  ssmd.hit <- log.ssmd.hit <- rep(0,times);
  result <- matrix(0,nrow=length(Decrease.fold),ncol=6)

  for(i in 1:length(Decrease.fold))
  {
    P.signal <- P.signal.o/Decrease.fold[i];
    for(j in 1:times)
    {
      N.random.spl <- sample(N.signal, sampling.num);
      P.random.spl <- sample(P.signal.o, sampling.num);
      N.random <- N.random.spl;
      P.random <- P.random.spl;
      ssmd.hit[j] <- gamma((sampling.num-1)/2)/gamma((sampling.num-2)/2)*sqrt(2/(sampling.num-1))*mean(P.random-median(N.random))/sd(P.random-median(N.random));


      N.random <- log2(N.random.spl);
      P.random <- log2(P.random.spl);
      log.ssmd.hit[j] <- gamma((sampling.num-1)/2)/gamma((sampling.num-2)/2)*sqrt(2/(sampling.num-1))*mean(P.random-median(N.random))/sd(P.random-median(N.random));
    }
    result[i,] <- c(mean(ssmd.hit),quantile((ssmd.hit),0.025),quantile((ssmd.hit),0.975),mean(log.ssmd.hit),quantile((log.ssmd.hit),0.025),quantile((log.ssmd.hit),0.975))
  }

  result <- cbind(Decrease.ratio,rep(sampling.num,length(Decrease.fold)),result)
  colnames(result) <- c("effect size","sampling number","mean_ssmd","95%CI_lower","95%CI_higher","mean_log2ssmd","95%CI_log2lower","95%CI_log2higher")
  if(is.null(dic.output)) dic.output <- getwd()
  write.csv(result, file=paste0(dic.output,"/ssmd_est.csv"), row.names=FALSE);
}


#####################################################
# REALTIME ANALYSIS-ONE DRUG FILE
#####################################################
#' @title Compound Analysis
#' @description This function creates a .png file that contains one Boxplot, SSMD, and Heatmap detailing the results of a drug compound plate.
#'
#' @param Pos.file A .xml file derived from your positive control condition, i.e., that with the highest value.
#' @param Neg.file A .xml file derived from your negative control condition, i.e., that with the lowest value.
#' @param Drug.file A .xml file derived from single drug treatment condition for comparison to controls.
#' @param well.label.file A .csv denoting the identity of each well of a 96-well plate, by condition or concentration.
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @param bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @param dic.output Output directory
#' @return  A .png file including: Boxplot, SSMD, and Heatmap detailing the drug compound plate.
#' @import grDevices
# @importFrom  pracma polyfit polyval
#' @import quadprog
#' @import pracma
#' @import graphics
#' @import stats
#' @import utils
#' @export

REAL.TIME <- function(Pos.file=NULL, Neg.file=NULL, Drug.file=NULL, well.label.file=NULL, Roy.file=NULL,bg.ratio=NULL,dic.output=NULL)
{
  if(is.null(well.label.file)) well.label.file <- system.file("extdata","Well.label.csv",package = "ARQiv")
  if(is.null(Pos.file)) Pos.file <- system.file("extdata","Positive.Control.xml",package = "ARQiv")
  if(is.null(Neg.file)) Neg.file <- system.file("extdata","Negative.Control.xml",package = "ARQiv")
  if(is.null(Drug.file)) Drug.file <- system.file("extdata","Drug.xml",package = "ARQiv")
  if(!is.null(Roy.file)){
    bg.ratio <- SIGNAL(Roy.file)
  }else{
    if(is.null(bg.ratio)) bg.ratio <- 3946
  }

  well.label1 <- read.csv(well.label.file,header = FALSE);

  neg_vector <- as.matrix(well.label1[3:10,2:13])
  neg_vector <- as.vector(t(neg_vector))

  pos_vector <- as.matrix(well.label1[15:22,2:13])
  pos_vector <- as.vector(t(pos_vector))

  compound_vector <- as.matrix(well.label1[27:34,2:13])
  compound_vector <- as.vector(t(compound_vector))

  well.label2 <- data.frame(Nplate=neg_vector,Pplate=pos_vector,Concentration=compound_vector)

  well.label <- well.label2

  well.name <- c(paste0("A", 1:12), paste0("B", 1:12), paste0("C", 1:12), paste0("D", 1:12), paste0("E", 1:12), paste0("F", 1:12), paste0("G", 1:12), paste0("H", 1:12));


  # NEGATIVE CONTROL SIGNAL EXTRACTION
  Neg.XML <- ParseXML(xmlfile=Neg.file);
  Neg.signal <- as.matrix(Neg.XML$signal[ ,1:9]);
  Nmax.signal <- apply(Neg.signal, 1, max);
  Neg.signal[Neg.signal<=bg.ratio] <- 0;
  Nwell.signal <- apply(Neg.signal, 1, sum);
  Nwell.signal <- as.numeric(apply(cbind(Nwell.signal, Nmax.signal), 1, max));
  s.neg <- Nwell.signal;
  s.neg <- s.neg[well.label$Nplate=="-"];

  # POSITIVE CONTROL SIGNAL EXTRACTION
  Pos.XML <- ParseXML(xmlfile=Pos.file);
  Pos.signal <- as.matrix(Pos.XML$signal[ ,1:9]);
  Pmax.signal <- apply(Pos.signal, 1, max);
  Pos.signal[Pos.signal<=bg.ratio] <- 0;
  Pwell.signal <- apply(Pos.signal, 1, sum);
  Pwell.signal <- as.numeric(apply(cbind(Pwell.signal, Pmax.signal), 1, max));
  s.pos <- Pwell.signal;
  ####plan2
  s200 <- s.pos[well.label$Pplate=="200"];
  s100 <- s.pos[well.label$Pplate=="100"];
  s50 <- s.pos[well.label$Pplate=="50"];
  s25 <- s.pos[well.label$Pplate=="25"];
  s12.5 <- s.pos[well.label$Pplate=="12.5"];
  s6.25 <- s.pos[well.label$Pplate=="6.25"];
  SSMD.s200.log <- (gamma(length(s200)/2-1/2)/gamma(length(s200)/2-1))*sqrt(2/(length(s200)-1))*mean(log(s200, base=2)-median(log(s.neg, base=2)))/sd(log(s200, base=2)-median(log(s.neg, base=2)));
  SSMD.s100.log <- (gamma(length(s100)/2-1/2)/gamma(length(s100)/2-1))*sqrt(2/(length(s100)-1))*mean(log(s100, base=2)-median(log(s.neg, base=2)))/sd(log(s100, base=2)-median(log(s.neg, base=2)));
  SSMD.s50.log <- (gamma(length(s50)/2-1/2)/gamma(length(s50)/2-1))*sqrt(2/(length(s50)-1))*mean(log(s50, base=2)-median(log(s.neg, base=2)))/sd(log(s50, base=2)-median(log(s.neg, base=2)));
  SSMD.s25.log <- (gamma(length(s25)/2-1/2)/gamma(length(s25)/2-1))*sqrt(2/(length(s25)-1))*mean(log(s25, base=2)-median(log(s.neg, base=2)))/sd(log(s25, base=2)-median(log(s.neg, base=2)));
  SSMD.s12.5.log <- (gamma(length(s12.5)/2-1/2)/gamma(length(s12.5)/2-1))*sqrt(2/(length(s12.5)-1))*mean(log(s12.5, base=2)-median(log(s.neg, base=2)))/sd(log(s12.5, base=2)-median(log(s.neg, base=2)));
  SSMD.s6.25.log <- (gamma(length(s6.25)/2-1/2)/gamma(length(s6.25)/2-1))*sqrt(2/(length(s6.25)-1))*mean(log(s6.25, base=2)-median(log(s.neg, base=2)))/sd(log(s6.25, base=2)-median(log(s.neg, base=2)));
  group_select <-which.max(c(SSMD.s200.log, SSMD.s100.log, SSMD.s50.log, SSMD.s25.log, SSMD.s12.5.log, SSMD.s6.25.log));
  s.con <- c("200","100","50","25","12.5","6.25")[group_select]
  s.pos <- s.pos[well.label$Pplate==s.con];
  ####

  # DRUG SIGNAL EXTRACTION
  Drug.XML <- ParseXML(xmlfile=Drug.file);
  Drug.signal <- as.matrix(Drug.XML$signal[ ,1:9]);
  Dmax.signal <- apply(Drug.signal, 1, max);
  Drug.signal[Drug.signal<=bg.ratio] <- 0;
  Dwell.signal <- apply(Drug.signal, 1, sum);
  Dwell.signal <- as.numeric(apply(cbind(Dwell.signal, Dmax.signal), 1, max));

  s4 <- Dwell.signal[well.label$Concentration=="4"];
  s2 <- Dwell.signal[well.label$Concentration=="2"];
  s1 <- Dwell.signal[well.label$Concentration=="1"];
  s0.5 <- Dwell.signal[well.label$Concentration=="0.5"];
  s0.25 <- Dwell.signal[well.label$Concentration=="0.25"];
  s0.125 <- Dwell.signal[well.label$Concentration=="0.125"];

  # GENERATION BOXPLOT
  if(is.null(dic.output)) dic.output <- getwd()
  png(file=paste0(dic.output,"/",Drug.XML$name,".plot.png"), width=1800, height=960);
  par(fig=c(0.02,0.31,0.05,1), font=2, lwd=3);
  boxplot.signal <- data.frame(rbind(cbind(s4/bg.ratio, rep("A", length(s4))), cbind(s2/bg.ratio, rep("B", length(s2))), cbind(s1/bg.ratio, rep("C", length(s1))), cbind(s0.5/bg.ratio, rep("D", length(s0.5))), cbind(s0.25/bg.ratio, rep("E", length(s0.25))), cbind(s0.125/bg.ratio, rep("F", length(s0.125))), cbind(s.neg/bg.ratio, rep("G", length(s.neg))), cbind(s.pos/bg.ratio, rep("H", length(s.pos)))));
  boxplot.signal <- data.frame(as.numeric(as.character(boxplot.signal[,1])), boxplot.signal[,2]);
  colnames(boxplot.signal) <- c("signal", "label");
  boxplot.par <- boxplot(signal ~ label, boxplot.signal, notch=TRUE, outline = FALSE, whisklty = 0, staplelty = 0, ylim=c(0,16), ylab="", names=rep("", 8), font.axis=2, font.lab=2, cex.axis=2, cex.lab=2, cex.main=2, lwd=3, yaxt='n', boxwex=0.5, plot=FALSE);
  ymax <- (max(boxplot.par$stats[4,]) %/% 8 +1) * 8;
  boxplot.par <- boxplot(signal ~ label, boxplot.signal, notch=TRUE, outline = FALSE, whisklty = 0, staplelty = 0, ylim=c(0,ymax), ylab="", names=rep("", 8), font.axis=2, font.lab=2, cex.axis=2, cex.lab=2, cex.main=2, lwd=3, yaxt='n', boxwex=0.5);
  axis(1, at=c(1:8), labels=c("4","2","1","0.5","0.25","0.125","(-)","(+)"), font=2, cex.axis = 2, lwd=3, las=3);
  axis(2, cex.axis = 2, lwd=3, font=2);
  mtext("YFP/Background", side=2, line=4, cex=2);
  mtext("Concentration", side=1, line=6, cex=2);
  s.median <- c(median(s4/bg.ratio), median(s2/bg.ratio), median(s1/bg.ratio), median(s0.5/bg.ratio), median(s0.25/bg.ratio), median(s0.125/bg.ratio));
  rank.polyfit <- polyfit(c(1:6), s.median, 2);
  rank.polyfit.value <- polyval(rank.polyfit, c(0:250)*0.02+1);
  lines(c(0:250)*0.02+1, rank.polyfit.value, lwd=3);

  # SSMD SCORE (HIT SELECTION) CALCULATION
  SSMD.pos.log <- (gamma(length(s.pos)/2-1/2)/gamma(length(s.pos)/2-1))*sqrt(2/(length(s.pos)-1))*mean(log(s.pos, base=2)-median(log(s.neg, base=2)))/sd(log(s.pos, base=2)-median(log(s.neg, base=2)));
  SSMD.neg.log <- (gamma(length(s.neg)/2-1/2)/gamma(length(s.neg)/2-1))*sqrt(2/(length(s.neg)-1))*mean(log(s.neg, base=2)-median(log(s.neg, base=2)))/sd(log(s.neg, base=2)-median(log(s.neg, base=2)));
  SSMD.s4.log <- (gamma(length(s4)/2-1/2)/gamma(length(s4)/2-1))*sqrt(2/(length(s4)-1))*mean(log(s4, base=2)-median(log(s.neg, base=2)))/sd(log(s4, base=2)-median(log(s.neg, base=2)));
  SSMD.s2.log <- (gamma(length(s2)/2-1/2)/gamma(length(s2)/2-1))*sqrt(2/(length(s2)-1))*mean(log(s2, base=2)-median(log(s.neg, base=2)))/sd(log(s2, base=2)-median(log(s.neg, base=2)));
  SSMD.s1.log <- (gamma(length(s1)/2-1/2)/gamma(length(s1)/2-1))*sqrt(2/(length(s1)-1))*mean(log(s1, base=2)-median(log(s.neg, base=2)))/sd(log(s1, base=2)-median(log(s.neg, base=2)));
  SSMD.s0.5.log <- (gamma(length(s0.5)/2-1/2)/gamma(length(s0.5)/2-1))*sqrt(2/(length(s0.5)-1))*mean(log(s0.5, base=2)-median(log(s.neg, base=2)))/sd(log(s0.5, base=2)-median(log(s.neg, base=2)));
  SSMD.s0.25.log <- (gamma(length(s0.25)/2-1/2)/gamma(length(s0.25)/2-1))*sqrt(2/(length(s0.25)-1))*mean(log(s0.25, base=2)-median(log(s.neg, base=2)))/sd(log(s0.25, base=2)-median(log(s.neg, base=2)));
  SSMD.s0.125.log <- (gamma(length(s0.125)/2-1/2)/gamma(length(s0.125)/2-1))*sqrt(2/(length(s0.125)-1))*mean(log(s0.125, base=2)-median(log(s.neg, base=2)))/sd(log(s0.125, base=2)-median(log(s.neg, base=2)));

  SSMD.log <- c(SSMD.s4.log, SSMD.s2.log, SSMD.s1.log, SSMD.s0.5.log, SSMD.s0.25.log, SSMD.s0.125.log, SSMD.neg.log, SSMD.pos.log);

  # SSMD SCORE PLOT
  par(fig=c(0.33,0.64,0.05,1), font=2, lwd=3, new=TRUE)
  plot(SSMD.log, ylim=c(-2,3), xlim=c(0.5, 8.5), pch='*', cex=4, xaxt='n', yaxt='n', ylab='', xlab='');
  axis(1, at=c(1:8), labels=c("4","2","1","0.5","0.25","0.125","(-)", "(+)"), font=2, cex.axis = 2, lwd=3, las=3);
  axis(2, at=c(-4:6)*0.5, labels=c(-4:6)*0.5 , cex.axis = 2, lwd=3, font=2);
  text(8,(SSMD.log[8]+0.2), c("200","100","50","25","12.5","6.25")[group_select],col="red",cex=1)
  abline(h=0, col="black", lwd=3, lty=3);
  abline(h=c(-0.25, 0.25), col="blue", lwd=3, lty=3);
  abline(h=c(-0.5, 0.5), col="gold", lwd=3, lty=3);
  abline(h=c(-1, 1), col="pink", lwd=3, lty=3);
  abline(h=c(-1.645, 1.645), col="green", lwd=3, lty=3);
  abline(h=c(-2, 2), col="red", lwd=3, lty=3);
  mtext("SSMD Score", side=2, line=4, cex=2);
  mtext("Concentration", side=1, line=6, cex=2);

  # GENERATE HEATMAP
  s.colorcode <- Dwell.signal[1:96]/(median(Dwell.signal[1:96])+1.5*IQR(Dwell.signal[1:96]));
  s.colorcode[s.colorcode>1] <- 1;
  s.x <- 1:96 %% 12;
  s.x[s.x==0] <- 12;
  s.y <- ceiling(1:96/12);

  # COLORRAMP PRODUCES CUSTOM PALETTES, BUT NEEDS VALUES BETWEEN 0 AND 1
  colorFunction <- colorRamp(c("black", "tan4", "darkgoldenrod"));

  # APPLY COLORRAMP AND SWITCH TO HEXADECIMAL REPRESENTATION
  s.colorcode.matrix <- colorFunction(s.colorcode);
  sColors <- rgb(s.colorcode.matrix, maxColorValue=255);

  # LET'S PLOT
  par(fig=c(0.61,0.99,0.3,1), font=1, lwd=2, new=TRUE);
  plot(x=s.x*1.5, y=(9-s.y)*0.94, col=sColors, pch=19, cex=8.3, xlim=c(0.5,18), ylim=c(0.5,10.5), bty="n", xaxt='n', yaxt="n", xlab='', ylab='');

  # PLOT AUTO FLUORESCENCE WELL
  auto.Rnum <- 5;
  Drug.autoflrc <- array(0, dim=96);
  Drug.autoflrc[apply(Drug.signal > bg.ratio, 1, sum)>=auto.Rnum] <- 1;
  autoflrc.x <- s.x[Drug.autoflrc==1];
  autoflrc.y <- s.y[Drug.autoflrc==1];
  if (length(autoflrc.x)>0)
  {
    text(x=autoflrc.x*1.5, y=(9-autoflrc.y)*0.94, labels="A", col="red", cex=3, family="Courier New");
  }
  dev.off();
}

#####################################################
# REALTIME ANALYSIS-ONE DRUG FILE
#####################################################
#' @title Compound Analysis-Multiple
#' @description This function creates multiple .png files that each contain one Boxplot, SSMD, and Heatmap detailing the results of a specific drug compound plate.
#'
#' @param Pos.file A .xml file derived from your positive control condition, i.e., that with the highest value. This 'bookends' sets of 10 drug plates.
#' @param Neg.file A .xml file derived from your negative control condition, i.e., that with the lowest value. This 'bookends' sets of 10 drug plates.
#' @param multi.Drug.file Multiple .xml files of drug treatment plates for for comparison to bordering controls.
#' @param well.label.file A .csv denoting the identity of each well of a 96-well plate, by condition or concentration.
#' @param Roy.file A .xml file of raw signal values from a 96-well plate of non-transgenic Roy larvae.
#' @param bg.ratio Numerical minimal 'signal' cutoff value (arbitrary units) applied to isolate fluorescent/luminescent signal above this value.
#' @param dic.output Output directory
#' @return multiple .png files which each include: Boxplot, SSMD, and Heatmap detailing a specific drug compound plate.
#' @import grDevices
#' @importFrom  pracma polyfit polyval
#' @import graphics
#' @import stats
#' @import utils
MULTI.REAL.TIME <- function(Pos.file=NULL, Neg.file=NULL, multi.Drug.file=NULL, well.label.file=NULL, Roy.file=NULL,bg.ratio=NULL,dic.output=NULL)
{
  for(i in 1:length(multi.Drug.file)){
    REAL.TIME(Pos.file =Pos.file, Roy.file=Roy.file, Neg.file=Neg.file, Drug.file=multi.Drug.file[i], well.label.file=Well.label.file,dic.output=dic.output)

  }
}

#####################################################
# GUI WINDOWS
#####################################################
#' @title GUI Window
#' @description The function produces a GUI window for users' convenience
#'
#' @return Boxplot SSMD and heatmap
#' @import grDevices
#' @import graphics
#' @import stats
#' @import utils
# @import gWidgets2
# @import RGtk2
# @import gWidgets2RGtk2
#' @export

GUI <- function(){
  options("guiToolkit"="RGtk2")
  win = gwindow("ARQiv HTS ANALYSIS", width=600, height=350)

  tbl = glayout(cont=win)
  #load files
  tbl[1,3:6, anchor=c(0,0)] <- glabel("PRE-SCREEN ASSAY OPTIMIZATION")
  tbl[2,1:8] <- gseparator(cont=tbl,horizontal = TRUE)
  tbl[3,2] <- glabel("INPUT")
  tbl[4,1] <- glabel("BACKGROUND")
  back_style <- gradio(c("FILE","VALUE"),horizontal = FALSE,cont=tbl)
  tbl[4:5,2] <- back_style
  tbl[4,3] <- gbutton(text = "bkgd", cont=tbl, handler = function(h,...){
    Roy.file <<- gfile()
    #  svalue(RoyFilename)<- ParseXML(Roy.file)$name
    svalue(RoyFilename)<- basename(Roy.file) })
  RoyFilename <- gedit(text = "", width=10, height=10, cont=tbl)
  tbl[4,5, anchor=c(0,0)] <- RoyFilename
  back_signal <- gspinbutton(from=0,to=6000,by=0.01, value=NULL, cont=tbl)
  tbl[5,3, anchor=c(0,0)] <- back_signal
  ############control
  tbl[6,1] <- glabel("CONTROL FILE")
  control_style <- gradio(c(".CSV",".XML"),selected=2,horizontal = FALSE,cont=tbl)
  tbl[6:7,2] <- control_style
  tbl[6,3] <- gbutton(text = "-/+", cont=tbl, handler = function(h,...){
    Neg.Pos.file <<- gfile()
    svalue(NegPosFilename)<- basename(Neg.Pos.file ) })
  NegPosFilename <- gedit(text = "", width=10, height=10, cont=tbl)
  tbl[6,5, anchor=c(0,0)] <- NegPosFilename
  tbl[7,3] <- gbutton(text = "+", cont=tbl, handler = function(h,...){
    Pos.file <<- gfile()
    # svalue(PosFilename)<- ParseXML(Pos.file)$name })
    svalue(PosFilename)<- basename(Pos.file)})
  PosFilename <- gedit(text = "", width=10, height=10, cont=tbl)
  tbl[7,5, anchor=c(0,0)] <- PosFilename
  tbl[8,3] <- gbutton(text = "-", cont=tbl, handler = function(h,...){
    Neg.file <<- gfile()
    #svalue(NegFilename)<- ParseXML(Neg.file)$name })
    svalue(NegFilename)<- basename(Neg.file)})
  NegFilename <- gedit(text = "", width=10, height=10, cont=tbl)
  tbl[8,5, anchor=c(0,0)] <- NegFilename

  ######effect size
  tbl[9,1] <- glabel("EFFECT SIZE")
  effect_style <- gradio(c("DEFAULT","CUSTOM"),selected=1,horizontal = FALSE,cont=tbl)
  tbl[9:10,2] <- effect_style
  effect_input <- gspinbutton(from=0,to=5,by=0.01, value=NULL, cont=tbl)
  tbl[9,3] <-glabel("(0.25,0.5,0.75,1)")
  tbl[10,3, anchor=c(0,0)] <- effect_input
  #### random sampling
  tbl[11,1] <- glabel("SAMPLING")
  tbl[11,2, anchor=c(0,0)] <- "SAMPLE#"
  sampling_number <- gspinbutton(from=0,to=96,by=1, value=0, cont=tbl)
  tbl[11,3, anchor=c(0,0)] <- sampling_number
  tbl[12,2, anchor=c(0,0)] <- "ITERATIONS"
  sampling_time <- gspinbutton(from=0,to=20000,by=1000, value=0, cont=tbl)
  tbl[12,3, anchor=c(0,0)] <- sampling_time
  tbl[2:13,6] <- gseparator(cont=tbl,horizontal = FALSE)


  ###OUTPUT
  tbl[3,7] <- glabel("OUTPUT")
  tbl[6,7] <- gbutton(text = "ASSAY DIRECTORY", cont=tbl, handler = function(h,...){
    dic.output <<- gfile(type="selectdir")
    svalue(DirectoryName)<- dic.output })
  DirectoryName <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[6,8, anchor=c(0,0)] <- DirectoryName
  ###
  tbl[4,7, anchor=c(0,0)] <- gbutton(text = "BACKGROUND", cont=tbl, handler = function(h,...){
    if(svalue(back_style)=="FILE"){
      svalue(bg.ratio)<- round(SIGNAL(Roy.file=Roy.file),2)
    }else svalue(bg.ratio) <- svalue(back_signal)
  })
  bg.ratio <-gedit(text = "", width=5, height=10, cont=tbl)
  tbl[4,8, anchor=c(0,0)] <- bg.ratio
  ###
  tbl[9,7, anchor=c(0,0)] <- gbutton(text = "SAMPLE SIZE", cont=tbl, handler = function(h,...){
    if(svalue(control_style)=="USE.XML"){
      if(svalue(effect_style)=="DEFAULT"){
        if(svalue(back_style)=="FILE"){
          SAMPLE(Pos.file =Pos.file, Roy.file=Roy.file,  Neg.file=Neg.file, Excel.file=NULL,effect.size =NULL, dic.output=dic.output)
        }else SAMPLE(Pos.file =Pos.file, Roy.file=NULL,bg.ratio=svalue(back_signal), Excel.file=NULL,  Neg.file=Neg.file, effect.size =NULL, dic.output=dic.output)

      }else{
        if(svalue(back_style)=="FILE"){
          SAMPLE(Pos.file =Pos.file, Roy.file=Roy.file,  Neg.file=Neg.file, Excel.file=NULL, effect.size =svalue(effect_input), dic.output=dic.output)
        }else SAMPLE(Pos.file =Pos.file, Roy.file=NULL,bg.ratio==svalue(back_signal),  Neg.file=Neg.file, Excel.file=NULL, effect.size =svalue(effect_input), dic.output=dic.output)
      }
    }else{
      if(svalue(effect_style)=="DEFAULT"){
        if(svalue(back_style)=="FILE"){
          SAMPLE(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=Roy.file, effect.size =NULL, dic.output=dic.output)
        }else  SAMPLE(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file,Roy.file=NULL,bg.ratio=svalue(back_signal),  effect.size =NULL, dic.output=dic.output)
      }else{
        if(svalue(back_style)=="FILE"){
          SAMPLE(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=Roy.file,  effect.size =svalue(effect_input), dic.output=dic.output)
        }else SAMPLE(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=NULL,bg.ratio==svalue(back_signal), effect.size =svalue(effect_input), dic.output=dic.output)
      }
    }
  })
  tbl[9,8, anchor=c(0,0)] <- glabel("(.CSV)")
  ###
  tbl[7,7, anchor=c(0,0)] <- gbutton(text = "Z'-FACTOR", cont=tbl, handler = function(h,...){
    #  Roy.file <<- gfile()
    if(svalue(control_style)==".XML"){
      if(svalue(back_style)=="FILE"){
        svalue(zfactor_1)<- ZFACTOR(Pos.file =Pos.file, Neg.file=Neg.file,Excel.file=NULL,Roy.file=Roy.file, bg.ratio=NULL , dic.output=dic.output)
      }else svalue(zfactor_1)<- ZFACTOR(Pos.file =Pos.file, Neg.file=Neg.file,Excel.file=NULL,Roy.file=NULL, bg.ratio=svalue(back_signal) , dic.output=dic.output)
    }else{
      if(svalue(back_style)=="FILE"){
        svalue(zfactor_1)<- ZFACTOR(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file,Roy.file=Roy.file, bg.ratio=NULL , dic.output=dic.output)
      }else svalue(zfactor_1)<- ZFACTOR(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file,Roy.file=NULL, bg.ratio=svalue(back_signal) , dic.output=dic.output)

    }
  })
  zfactor_1 <-gedit(text = "", width=15, height=10, cont=tbl)
  tbl[7,8, anchor=c(0,0)] <-zfactor_1
  ####
  tbl[8,7, anchor=c(0,0)] <- gbutton(text = "SSMD QC", cont=tbl, handler = function(h,...){
    #  Roy.file <<- gfile()
    if(svalue(control_style)==".XML"){
      if(svalue(back_style)=="FILE"){
        svalue(ssmd.qc_1)<- SSMD.QC(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL, Roy.file=Roy.file, bg.ratio=NULL ,dic.output=dic.output)
      }else svalue(ssmd.qc_1)<- SSMD.QC(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL, Roy.file=NULL, bg.ratio=svalue(back_signal) ,dic.output=dic.output)
    }else{
      if(svalue(back_style)=="FILE"){
        svalue(ssmd.qc_1)<- SSMD.QC(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=Roy.file, bg.ratio=NULL ,dic.output=dic.output)
      }else svalue(ssmd.qc_1)<- SSMD.QC(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=NULL, bg.ratio=svalue(back_signal) ,dic.output=dic.output)

    }
  })
  ssmd.qc_1 <-gedit(text = "", width=15, height=10, cont=tbl)
  tbl[8,8, anchor=c(0,0)] <-ssmd.qc_1
  ####
  tbl[10,7, anchor=c(0,0)] <- gbutton(text = "SSMD ESTIMATE", cont=tbl, handler = function(h,...){
    #  Roy.file <<- gfile()
    if(svalue(sampling_number)!=0&svalue(sampling_time)!=0){
      if(svalue(back_style)=="FILE"){
        if(svalue(control_style)==".XML"){
          if(svalue(back_style)=="COSTOM"){
            SSMD.EST(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL, Roy.file=Roy.file, effect.size=svalue(effect_input), sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
          }else SSMD.EST(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL, Roy.file=Roy.file, effect.size=NULL, sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
        }else{
          if(svalue(back_style)=="FILE"){
            SSMD.EST(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=Roy.file, effect.size=svalue(effect_input), sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
          }else SSMD.EST(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=Roy.file, effect.size=NULL, sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
        }
      }else{
        if(svalue(control_style)==".XML"){
          if(svalue(back_style)=="COSTOM"){
            SSMD.EST(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL, Roy.file=NULL, bg.ratio=svalue(back_signal), effect.size=svalue(effect_input), sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
          }else SSMD.EST(Pos.file =Pos.file,Neg.file=Neg.file,Excel.file=NULL,Roy.file=NULL, bg.ratio=svalue(back_signal), effect.size=NULL, sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
        }else{
          if(svalue(back_style)=="FILE"){
            SSMD.EST(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=NULL, bg.ratio=svalue(back_signal), effect.size=svalue(effect_input), sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
          }else SSMD.EST(Pos.file =NULL,Neg.file=NULL,Excel.file =Neg.Pos.file, Roy.file=NULL, bg.ratio=svalue(back_signal), effect.size=NULL, sampling.num=svalue(sampling_number),times=svalue(sampling_time),dic.output=dic.output)
        }
      }
    }
  })
  tbl[10,8, anchor=c(0,0)] <- glabel("(.CSV)")
  tbl[13,1:8] <- gseparator(cont=tbl,horizontal = TRUE)
  tbl[14,1:8] <- gseparator(cont=tbl,horizontal = TRUE)
  tbl[15,1:8] <- gseparator(cont=tbl,horizontal = TRUE)
  ####compond analysis
  tbl[16,3:7, anchor=c(0,0)] <- glabel("COMPOUND ANALYSIS")
  tbl[17,1:8] <- gseparator(cont=tbl,horizontal = TRUE)
  tbl[18,2] <- glabel("INPUT")
  tbl[19,1] <- glabel("BACKGROUND")
  back_style1 <- gradio(c("FILE","VALUE"),horizontal = FALSE,cont=tbl)
  tbl[19:20,2] <- back_style1
  tbl[19,3] <- gbutton(text = "bkgd", cont=tbl, handler = function(h,...){
    Roy.file1 <<- gfile()
    svalue(RoyFilename1)<- ParseXML(Roy.file1)$name })
  RoyFilename1 <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[19,5, anchor=c(0,0)] <- RoyFilename1
  back_signal1 <- gspinbutton(from=0,to=6000,by=0.01, value=NULL, cont=tbl)
  tbl[20,3, anchor=c(0,0)] <- back_signal1
  ####control
  tbl[21,1] <- glabel("CONTROL.XML")
  tbl[21,2] <- gbutton(text = "+", cont=tbl, handler = function(h,...){
    Pos.file1 <<- gfile()
    svalue(PosFilename1)<- basename(Pos.file1) })
  PosFilename1 <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[21,3, anchor=c(0,0)] <- PosFilename1
  tbl[23,2] <- gbutton(text = "-", cont=tbl, handler = function(h,...){
    Neg.file1 <<- gfile()
    svalue(NegFilename1)<- basename(Neg.file1) })
  NegFilename1 <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[23,3, anchor=c(0,0)] <- NegFilename1
  ###well
  tbl[24,1] <- glabel("WELL LABEL FILE")
  tbl[24,2] <- gbutton(text = ".CSV", cont=tbl, handler = function(h,...){
    Well.label.file <<- gfile()
    svalue(WellLabelFilename)<- basename(Well.label.file) })
  WellLabelFilename <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[24,3, anchor=c(0,0)] <- WellLabelFilename
  ####drug files
  tbl[19,7] <- gbutton(text = "DRUG FILES(.XML)", cont=tbl, handler = function(h,...){
    r<- gfile(multi=T)
    r<- unlist(strsplit(r, "\n"))
    multi.Drug.file <<- r
    for(i in 1:length(r)) {
     # r[i] <- ParseXML(r[i])$name
      r[i] <- basename(r[i])
    }
    svalue(DrugXMLFilenames)<- paste(r, sep="\n")
    focus(DrugXMLFilenames)<- T })
  DrugXMLFilenames <- gtext(text = "", width=15, height=110, cont=tbl)
  tbl[19:21,8, anchor=c(0,0)] <- DrugXMLFilenames
  tbl[22,4:8] <- gseparator(cont=tbl)
  tbl[22:25,4] <- gseparator(cont=tbl,horizontal = FALSE)
  tbl[23,5] <- glabel("OUTPUT")
  tbl[24,5] <- gbutton(text = "COMPOUND DIRECTORY", cont=tbl, handler = function(h,...){
    dic.output1 <<- gfile(type="selectdir")
    svalue(DirectoryName1)<- dic.output1 })
  DirectoryName1 <- gedit(text = "", width=15, height=10, cont=tbl)
  tbl[24,6:7, anchor=c(0,0)] <- DirectoryName1
  #ANALYSIS
  tbl[25,5:6, anchor=c(0,0)] <- gbutton(text = "ANALYSIS", cont=tbl, handler = function(h,...){
    if(svalue(back_style1)=="FILE"){
      MULTI.REAL.TIME(Pos.file =Pos.file1, Roy.file=Roy.file1,bg.ratio = NULL, Neg.file=Neg.file1, multi.Drug.file=multi.Drug.file, well.label.file=Well.label.file,dic.output=dic.output1)
    }else MULTI.REAL.TIME(Pos.file =Pos.file1, Roy.file=NULL,bg.ratio = svalue(back_signal1), Neg.file=Neg.file1, multi.Drug.file=multi.Drug.file, well.label.file=Well.label.file,dic.output=dic.output1)
  })

  tbl[26,1:9] <- gseparator(cont=tbl)
}

