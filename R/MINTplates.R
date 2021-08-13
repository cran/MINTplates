#' Encode a sequence
#' 
#' Encode the sequence into its corresponding license plate with given prefix (if given one). 
#' 
#' @param sequence The sequence being encoded.
#' @param prefix The prefix to use for the license plate.
#' @return The license plate it encodes to.
#' @export
#' @references 
#' Pliatsika V, Loher P, Telonis AG, Rigoutsos I. Bioinformatics (2016) <doi: 10.1093/bioinformatics/btw194>
#' @examples
#' enseq('AACCGGGCAGAAGCACCA','tRF')
enseq=function(sequence,prefix='')
{
  sequence=gsub('U','T',toupper(sequence))
  sequence=gsub("^\\s+|\\s+$", "",sequence)
  sequence=gsub("[\r\n\t ]", "",sequence)
  isseq=grep("^[ATCG]+$",sequence)
  if(length(isseq)==0)
    return('It is not a DNA/RNA sequence')
  encoded=paste(nchar(sequence),'-',sep = '')
  bk=FALSE
  repeat {
    s=sequence
    if(nchar(sequence)>5)
    {
      s=substr(s,1,5)
      sequence=substr(sequence,6,nchar(sequence))
    }
    else
    {
      bk=TRUE
    }
    encoded=paste(encoded,codes[codes$seq==s,'code'],sep = '')
    
    if(bk)
      break
  }
  if(prefix=='')
    return(encoded)
  else
    return(paste(prefix,encoded,sep = '-'))
}

#' Encode multiple sequences
#' 
#' Encode the sequences into their corresponding license plates with given prefix (if given one). 
#' 
#' @param sequences The sequences being encoded.
#' @param prefix The prefix to use for the license plate.
#' @return The license plates they encode to.
#' @export
#' @references 
#' Pliatsika V, Loher P, Telonis AG, Rigoutsos I. Bioinformatics (2016) <doi: 10.1093/bioinformatics/btw194>
#' @examples
#' seqs=c('AACCGGGCAGAAGCACCA','GAGCCCCAGTGGAACCACCA')
#' enseqs(seqs,'tRF')
enseqs=function(sequences,prefix='')
{
  plates=list()
  for(seq in sequences)
  {
    p=enseq(seq,prefix)
    plates=c(plates,list(p))
  }
  return(unlist(plates))
}

#' Decode a license-plate
#' 
#' Decode the license-plate using the lookup table. 
#' 
#' @param plate The license plate being decoded.
#' @return The sequence it decodes to.
#' @export
#' @references 
#' Pliatsika V, Loher P, Telonis AG, Rigoutsos I. Bioinformatics (2016) <doi: 10.1093/bioinformatics/btw194>
#' @examples
#' deseq('tRF-18-BS6PDFD2')
deseq=function(plate)
{
  plate=gsub("^\\s+|\\s+$", "",plate)
  plateLst=as.character(unlist(strsplit(plate, split = "-")))
  if(length(plateLst)==3)
    plateLst=plateLst[-1]
  
  isNum=grep("^[1234567890]+$",plateLst[1])
  isCode=grep("^[1234567890BDEFHIJKLMNOPQRSUVWXYZ]+$",plateLst[2])
  if(length(plateLst)!=2 | length(isNum)==0 | length(isCode)==0)
  {
    return('Bad format for plate')
  }
  
  
  Len=as.numeric(plateLst[1])
  code=plateLst[2]
  
  lastCodeLen=Len%%5

  CodeLen=5
  decoded=''
  
  bk=FALSE
  repeat {
    s=code
    if(nchar(code)>2)
    {
      s=substr(s,1,2)
      code=substr(code,3,nchar(code))
    }
    else
    {
      if(lastCodeLen!=0)
         CodeLen=lastCodeLen
      bk=TRUE
    }

    index=which(codes$length==CodeLen & codes$code==s)
    decoded=paste(decoded,codes[index,'seq'],sep = '')
    if(bk)
      break
  }
  return(decoded)
}

#' Decode multiple license-plates
#' 
#' Decode the license-plates using the lookup table. 
#' 
#' @param plates The license plates being decoded.
#' @return The sequences they decodes to.
#' @export
#' @references 
#' Pliatsika V, Loher P, Telonis AG, Rigoutsos I. Bioinformatics (2016) <doi: 10.1093/bioinformatics/btw194>
#' @examples
#' plates=c('tRF-18-BS6PDFD2','tRF-20-51K36D26')
#' deseqs(plates)
deseqs=function(plates)
{
  seqs=list()
  for(p in plates)
  {
    p=deseq(p)
    seqs=c(seqs,list(p))
  }
  return(unlist(seqs))
}
