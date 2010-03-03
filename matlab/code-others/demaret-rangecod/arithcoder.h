#ifndef _ARITHCODER_H_
#define _ARITHCODER_H_

#include "qsmodel.h"
#include "rangecod.h"
#include "bitstr.h"

//#define BASE_2_LOG_FREQ_COUNT   12 // base2 log of total frequency count
#define BASE_2_LOG_FREQ_COUNT   13 // base2 log of total frequency count
#define RESCALING_INTERVAL      2000


class ArithmCoder
{
 public:
  // Constructor
  ArithmCoder(unsigned int NbOfSymbols);

  // Destructor
  virtual ~ArithmCoder();

  // Accessors
  unsigned int GetNbEncodedSymb();

  // Encoding functions
  void StartEncoding();
  void EncodeValue(unsigned int value, WBitStream& OutBitStream);
  void EncodeValueContext(unsigned int ContextMax,unsigned int value, WBitStream& OutBitStream);
  void EncodeValueContext(unsigned int ContextMin, unsigned int LocalMax,unsigned int value, WBitStream& OutBitStream);
  void EncodeValueCombination(double proba0,unsigned int value, WBitStream& OutBitStream);
  void EncodeValueCombination(double cum_proba0, double cum_proba1, unsigned int value,
                                                 WBitStream& OutBitStream);
 void DoneEncoding(WBitStream& OutBitStream);


  // Decoding functions
  void StartDecoding(RBitStream& ro_InBitStream);
  void StartDecodingDeterministic(RBitStream& InBitStream);

  bool DecodeValue(unsigned int& rui_Value, RBitStream& InBitStream);
  bool DecodeValue(int* Value, RBitStream& InBitStream);
  bool DecodeValueContext(unsigned int ContextMax,int* value, RBitStream& InBitStream);
  bool DecodeValueContext(unsigned int ContextMin,unsigned int ContextMax,int* value, RBitStream& InBitStream);
  bool DecodeValueCombination(double proba0, int* Value, RBitStream& InBitStream);
  bool DecodeCombinationMulti(double* proba, int* Value, RBitStream& InBitStream);


  void DoneDecoding(RBitStream& ro_InBitStream);


 protected:
  unsigned int NbSymbols;  // Does not include EOF
  unsigned int NbEncodedSymbols;
  // Data
  rangecoder rc;
  qsmodel qsm;
  int syfreq, ltfreq;

};

#endif
