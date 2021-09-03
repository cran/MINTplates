#' Encode sequences
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
  plates=c()
  for(sequence in sequences)
  {
    p='It is not a DNA/RNA sequence'
    sequence=gsub('U','T',toupper(sequence))
    sequence=gsub("^\\s+|\\s+$", "",sequence)
    sequence=gsub("[\r\n\t ]", "",sequence)
    isseq=grep("^[ATCG]+$",sequence)
    if(length(isseq)==0)
    {
      plates=c(plates,p)
      next
    }
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
      encoded=paste(encoded,dat$codes[dat$codes$seq==s,'code'],sep = '')
      
      if(bk)
        break
    }
    if(prefix=='')
      p=(encoded)
    else
      p=(paste(prefix,encoded,sep = '-'))
    plates=c(plates,p)
  }
  return(plates)
}


#' Decode license-plates
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
  seqs=c()
  for(plate in plates)
  {
    p='Bad format for plate'
    
    
    plate=gsub("^\\s+|\\s+$", "",plate)
    plateLst=as.character(unlist(strsplit(plate, split = "-")))
    if(length(plateLst)==3)
      plateLst=plateLst[-1]
    
    isNum=grep("^[1234567890]+$",plateLst[1])
    isCode=grep("^[1234567890BDEFHIJKLMNOPQRSUVWXYZ]+$",plateLst[2])
    if(length(plateLst)!=2 | length(isNum)==0 | length(isCode)==0)
    {
      seqs=c(seqs,(p))
      next
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
      
      index=which(dat$codes$length==CodeLen & dat$codes$code==s)
      decoded=paste(decoded,dat$codes[index,'seq'],sep = '')
      if(bk)
        break
    }
    p=(decoded)
    
    seqs=c(seqs,(p))
  }
  return((seqs))
}

#' Obtain the type of tRFs
#' 
#' Obtain the type of tRFs with the tRF sequences. 
#' 
#' @param sequences tRF sequences.
#' @return The type of tRFs.
#' @export
#' @references 
#' Loher P, Telonis AG, Rigoutsos I. Sci Rep (2017) <doi: 10.1038/srep41184>
#' @examples
#' sequences=c('TCCCTGGTGGTCTAGTGGTTAGGATTCGGC','AAAAATTTTGGTGCAACTCCAAATAAAA')
#' tRFtype(sequences)
tRFtype=function(sequences)
{
  return(dat$LookupTable[sequences,'category'])
}

#' Determine whether the tRFs are exclusive to tRNA space
#' 
#' Determine whether the tRFs are exclusive to tRNA space with the tRF sequences. 
#' 
#' @param sequences tRF sequences.
#' @return Whether the tRFs are exclusive to tRNA space.
#' @export
#' @references 
#' Loher P, Telonis AG, Rigoutsos I. Sci Rep (2017) <doi: 10.1038/srep41184>
#' @examples
#' sequences=c('TCCCTGGTGGTCTAGTGGTTAGGATTCGGC','TCCCTGGTGGTCTAGTGGTTAGGATTCGGCG')
#' exclusive(sequences)
exclusive=function(sequences)
{
  return(dat$LookupTable[sequences,'exclusive'])
}

#' Obtain the tRNA source of a tRF
#' 
#' Obtain the tRNA source of a tRF with the tRF sequence. 
#' 
#' @param sequence tRF sequence.
#' @return Sources of the tRF.
#' @export
#' @references 
#' Loher P, Telonis AG, Rigoutsos I. Sci Rep (2017) <doi: 10.1038/srep41184>
#' @examples
#' sequence='TCCCTGGTGGTCTAGTGGTTAGGATTCGGC'
#' source_tRNA(sequence)
source_tRNA=function(sequence)
{
  trs=data.frame('amino_acid'=character(0),'anticodon'=character(0),'chr'=character(0))
  for (i in (1:nrow(dat$tRNAs))) {
    aa=dat$tRNAs[i,'amino_acid']
    seq=dat$tRNAs[i,'seq']
    tsSeq=substr(sequence,2,nchar(sequence))
    isIn=grep(tsSeq,seq)
    #print(isIn)
    if(length(isIn)!=0){
      b=data.frame(dat$tRNAs[i,c('amino_acid','anticodon','chr')])
      trs=rbind(trs,b)
    }
    
  }
  rownames(trs)=NULL
  return(trs)
}

#' Obtain the tRNA sources of tRFs
#' 
#' Obtain the tRNA sources of tRFs with the tRF sequences. 
#' 
#' @param sequences tRF sequences.
#' @return Sources of the tRFs.
#' @export
#' @references 
#' Loher P, Telonis AG, Rigoutsos I. Sci Rep (2017) <doi: 10.1038/srep41184>
#' @examples
#' sequences=c('TCCCTGGTGGTCTAGTGGTTAGGATTCGGC','TCCCTGGTGGTCTAGTGGCT','TCCCTGGTGGTCTAATGGTTA')
#' source_tRNA2(sequences)
source_tRNA2=function(sequences)
{
  tRNAs=list()
  for(sequence in sequences)
  {
    tRNA=''
    trs=source_tRNA(sequence)
    if(nrow(trs)==0)
    {
      tRNAs=c(tRNAs,list(NA))
      next
    }
    for(i in (1:nrow(trs)))
    {
      tRNA=paste(tRNA,trs[i,'amino_acid'],sep = '')
      tRNA=paste(tRNA,trs[i,'anticodon'],sep = '-')
      tRNA=paste(tRNA,trs[i,'chr'],sep = ' ')
      if(i<nrow(trs))
        tRNA=paste(tRNA,'; ',sep = '')
    }
    tRNAs=c(tRNAs,list(tRNA))
  }
  
  return(unlist(tRNAs))
}

#' Annotate a tRF sequence
#' 
#' Obtain tRF ID, type, whether exclusive to tRNA space, and tRNA sources of the tRF with its sequence. 
#' 
#' @param sequence tRF sequence.
#' @return tRF ID, type, whether exclusive to tRNA space, and tRNA sources of the tRF.
#' @export
#' @references 
#' Loher P, Telonis AG, Rigoutsos I. Sci Rep (2017) <doi: 10.1038/srep41184>
#' @examples
#' sequence='TCCCTGGTGGTCTAGTGGTTAGGATTCGGC'
#' annotate_tRF(sequence)
annotate_tRF=function(sequence)
{
  tRF=c('Sequence'=(sequence))
  tRF=c(tRF,list('tRF ID'=enseqs(sequence)))
  buf=dat$LookupTable[sequence,c('category','exclusive'),]
  tRF=c(tRF,list('Type'=buf[1,'category']))
  tRF=c(tRF,list('Exclusive to tRNA space'=buf[1,'exclusive']))
  tRF=c(tRF,list('Sources of tRNAs'=source_tRNA(sequence)))
  return(tRF)
}
