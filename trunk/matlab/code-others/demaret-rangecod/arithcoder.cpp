//  ARITHMETIC ENCODING ALGORITHM
#include "qsmodel.h"
#include "rangecod.h"
#include "bitstr.h"
#include "arithcoder.h"

using namespace std;


// Constructor with number of symbols
ArithmCoder::ArithmCoder(unsigned int NbSymb)
{
  NbSymbols = NbSymb;
  NbEncodedSymbols = 0;
}


// Destructor
ArithmCoder::~ArithmCoder()
{
}


unsigned int ArithmCoder::GetNbEncodedSymb()
{
  return NbEncodedSymbols;
}


void ArithmCoder::StartEncoding()
{
  // We need NbSymbols+1 because of EOF
  initqsmodel(&qsm,NbSymbols+1,BASE_2_LOG_FREQ_COUNT,RESCALING_INTERVAL,NULL,1);
  start_encoding(&rc,0);
}


void ArithmCoder::EncodeValue(unsigned int value, WBitStream& OutBitStream)
{
  qsgetfreq(&qsm,value,&syfreq,&ltfreq);
  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}


void ArithmCoder::EncodeValueContext(unsigned int LocalMax,unsigned int value, WBitStream& OutBitStream)
{
  //this line changes vs no context 	 
  qsgetfreqcont(&qsm,value,LocalMax,&syfreq,&ltfreq);

  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}




//same as previous : only values between localmin and localmax
void ArithmCoder::EncodeValueContext(unsigned int LocalMin, unsigned int LocalMax,unsigned int value, WBitStream& OutBitStream)
{
  //this line changes vs no context 	 
  qsgetfreqcont(&qsm,value,LocalMin,LocalMax,&syfreq,&ltfreq);

  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}


//This function encodes a combination of 0 and 1: 
// 0 with probability proba0 and 1 with probability (1-proba0)
void ArithmCoder::EncodeValueCombination(double proba0,unsigned int value, WBitStream& OutBitStream)
{
  //this line changes vs no context 	 
  qsgetfreqcomb(&qsm,value,proba0,&syfreq,&ltfreq);

  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}


//Chantier !!!!
//This function encode a given symbol according to a given probability
//cum_proba0 is the probability to be smaller than the symbol
//cum_proba1 is the probability to be bigger or equal than the symbol
void ArithmCoder::EncodeValueCombination(double cum_proba0, double cum_proba1, unsigned int value,
                                                 WBitStream& OutBitStream)
{
  //this line changes vs no context 	 
  qsgetfreqcombmulti(&qsm,cum_proba0,cum_proba1,&syfreq,&ltfreq);
  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  qsupdate(&qsm,value);
  NbEncodedSymbols++;
}




void ArithmCoder::DoneEncoding(WBitStream& OutBitStream)
{
  qsgetfreq(&qsm,NbSymbols,&syfreq,&ltfreq);
  encode_shift(&rc,syfreq,ltfreq,BASE_2_LOG_FREQ_COUNT,OutBitStream);
  done_encoding(&rc,OutBitStream);
  deleteqsmodel(&qsm);
  NbEncodedSymbols++;
}


void ArithmCoder::StartDecoding(RBitStream& InBitStream)
{
  // NbSymbols+1 because of EOF
  initqsmodel(&qsm,NbSymbols+1,BASE_2_LOG_FREQ_COUNT,RESCALING_INTERVAL,NULL,0);
  start_decoding(&rc,InBitStream);
}

void ArithmCoder::StartDecodingDeterministic(RBitStream& InBitStream)
{
  // only NbSymbols, because EOF does not need to be coded (the number of symbols is known)
  initqsmodel(&qsm,NbSymbols,BASE_2_LOG_FREQ_COUNT,RESCALING_INTERVAL,NULL,0);
  start_decoding(&rc,InBitStream);
}



bool ArithmCoder::DecodeValueContext(unsigned int LocalMax,int* Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  //Value = qsgetsymcont(&qsm,ltfreq,RestrictedMax);
  qsgetsymcont(&qsm,ltfreq,LocalMax,Value);
  int LocalValue = *Value;
  NbEncodedSymbols++;
  if (LocalValue == NbSymbols)  /* check for end-of-file */
    return false;

  qsgetfreqcont(&qsm,LocalValue,LocalMax,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,LocalValue);

  return true;
}


bool ArithmCoder::DecodeValueContext(unsigned int LocalMin, unsigned int LocalMax,int* Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  //Value = qsgetsymcont(&qsm,ltfreq,RestrictedMax);
  qsgetsymcont(&qsm,ltfreq,LocalMin,LocalMax,Value);
  int LocalValue = *Value;
  NbEncodedSymbols++;
  if (LocalValue == NbSymbols)  /* check for end-of-file */
    return false;

  qsgetfreqcont(&qsm,LocalValue,LocalMin,LocalMax,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,LocalValue);

  return true;
}


bool ArithmCoder::DecodeValueCombination(double proba0,int* Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  qsgetsymcomb(&qsm,ltfreq,proba0,Value);
  int LocalValue = *Value;
  NbEncodedSymbols++;
  qsgetfreqcomb(&qsm,LocalValue,proba0,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,LocalValue);

  return true;
}


//chantier:
//Remark : this works differently from the others, the probability is decoded,
//to get the value you have to have at your disposal the table of probabilities
bool ArithmCoder::DecodeCombinationMulti(double* proba, int* Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  int ltfreq_dec = ltfreq;
  qsgetcombmult(&qsm,ltfreq_dec,proba,Value,&syfreq,&ltfreq); // does also the job of qsgetfreqcomb
    
  int LocalValue = *Value;
  NbEncodedSymbols++;
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT); 
  qsupdate(&qsm,LocalValue); 

  return true;
}


bool ArithmCoder::DecodeValue(unsigned int& Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  Value = qsgetsym(&qsm,ltfreq);
  NbEncodedSymbols++;
  if (Value == NbSymbols)  /* check for end-of-file */
    return false;

  qsgetfreq(&qsm,Value,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,Value);

  return true;
}


bool ArithmCoder::DecodeValue(int* Value, RBitStream& InBitStream)
{
  ltfreq = decode_culshift(&rc,BASE_2_LOG_FREQ_COUNT,InBitStream);
  *Value = qsgetsym(&qsm,ltfreq);
  int LocalValue = *Value;
  NbEncodedSymbols++;
  if (*Value == NbSymbols)  /* check for end-of-file */
    return false;

  qsgetfreq(&qsm,LocalValue,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  qsupdate(&qsm,LocalValue);

  return true;
}


void ArithmCoder::DoneDecoding(RBitStream& InBitStream)
{
  qsgetfreq(&qsm,NbSymbols,&syfreq,&ltfreq);
  decode_update(&rc,syfreq,ltfreq,1<<BASE_2_LOG_FREQ_COUNT);
  done_decoding(&rc,InBitStream);
  deleteqsmodel(&qsm);
}
