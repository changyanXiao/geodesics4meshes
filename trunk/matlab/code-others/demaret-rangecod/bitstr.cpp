#include <iostream>
#include <fstream>

#include "bitstr.h"

using namespace std;

//Default constructor
BitStream::BitStream()
{
  Buffer = 0;
  tToStartBit();
  Counter = 0;
}


// Destructor
BitStream::~BitStream()
{
}


// Read the current bit in the stream
unsigned int BitStream::ReadCurrentBit()
{
  unsigned int bit = (Mask == (Mask & Buffer));
  return bit;
}


void BitStream::tToStartBit()
{
 // Current bit position initialized with the Most Significant Bit of the Buffer
  BitPos = sizeof(unsigned char)*8-1;

  // Mask placed on the Most Significant Bit of the Buffer
  Mask = 1 << BitPos;
}

// Position the current bit to the least significant bit in
// the buffer
void BitStream::tToEndBit()
{
  // Current bit position initialized with the Low Significant Bit of the Buffer
  BitPos = 0;
  // Mask placed on the Most Significant Bit of the Buffer
  Mask = 1;
}


// Write a bit in the stream
// Warning : Bit is uunsigned int but must be 0 or 1
void BitStream::WriteCurrentBit(unsigned int Bit)
{
  if (Bit == 0)
    Buffer &= (~Mask);
  else
    if (Bit == 1)
      Buffer |= Mask;
}


// Write a bit in the stream
void WriteCurrentBit(unsigned char Bit)
{
  WriteCurrentBit((unsigned int) Bit);
}


// Moving to the next buffer bit
unsigned int BitStream::ToNextBit()
{
  Counter++;
  if (BitPos == 0)
  {
    tToStartBit();
    return 1;
  }
  else
  {
    BitPos--;
    Mask >>= 1;
    return 0;
  }
}


// Moves to the previous buffer bit
unsigned int BitStream::ToPreviousBit()
{
  Counter--;
  if (BitPos == 7)
  {
    tToEndBit();
    return 1;
  }
  else
  {
    BitPos++;
    Mask <<= 1;
    return 0;
  }
}


// Returns counter of the bitstream
unsigned int BitStream::GetCounter()
{
  return Counter;
}


// Functions specific to WBitStream class




//Default constructor
WBitStream::WBitStream()
    : BitStream()
{
  pStream = 0;
}

// Purpose     : Constructor with the file name to open for writing
WBitStream::WBitStream(const char* pFileName)
  : BitStream()
{
  // Opends the stream (writing)
  pStream = new ofstream(pFileName,ios::binary);
  pStream->seekp(0);
}



// Destructor
WBitStream::~WBitStream()
{
  delete pStream;
}


// Save value to bitstream with ui_NbBitdsUsed bits (the size is checked)
// Return FALSE if an error occured
int WBitStream::SaveValue(unsigned int NbBitsUsed, unsigned int ValueToEncode, char* pc_ErrorMessage/*=0*/)
{
  int b_Error = (unsigned int)(1 << NbBitsUsed) <= ValueToEncode;
  if (b_Error)
  {
    cout <<  "Error in header saving";
    if (pc_ErrorMessage) cout << " :\n   While saving "<< *pc_ErrorMessage;
    cout << endl;
  }
  else
    WriteBits(ValueToEncode, NbBitsUsed);

  return !b_Error;
}


// Stuff the end of buffer with a given value
void WBitStream::Stuff(unsigned int Bit)
{
  do
  {
    WriteCurrentBit(Bit);
  }
  while (!ToNextBit());
  pStream->put(Buffer);
  Buffer = 0;
};



// Write a bit to the output stream
void WBitStream::WriteBit(unsigned int Bit)
{
  WriteCurrentBit(Bit);
  if (ToNextBit())
  {
    pStream->put(Buffer);
    Buffer = 0;
  }
}


// Write a bit to the ouput stream
void WBitStream::WriteBit(unsigned char Bit)
{
  WriteBit((unsigned int) Bit);
}



//      : Write bits on the file
// Warning we must have : NbBits <= sizeof(unsigned int)*8
void WBitStream::WriteBits(unsigned int Word, unsigned int NbBits)
{
  unsigned int WordMask = 1;
  unsigned int CurrentPos;
  unsigned int BitVal;

  for (CurrentPos = 0; CurrentPos < NbBits; CurrentPos++)
  {
    BitVal = ((Word & WordMask) ==  WordMask);
    WriteBit(BitVal);
    WordMask <<= 1;
  }
}


void WBitStream::WriteFloat(float Value)
{
  float* pValue = &Value;
  unsigned int *ptmp = (unsigned int *)pValue;
  WriteBits(*ptmp,sizeof(unsigned int)*8);
}


void WBitStream::WriteDouble(double d_Value)
{
  double* pValue = &d_Value;
  unsigned int *ptmp = (unsigned int *)pValue;
  WriteBits(*ptmp,sizeof(unsigned int)*8);
  ptmp++;
  WriteBits(*ptmp,sizeof(unsigned int)*8);
}


void WBitStream::WriteByte(unsigned char Byte)
{
  if (BitPos < 7)
  {
    Stuff(0);
    pStream->flush();
  }

  pStream->put(Byte);
  Counter += 8;
}


void WBitStream::WriteSignBits(int SignWord, unsigned int NbBits)
{
  if (SignWord >= 0) WriteBit((unsigned char)0);
  else WriteBit((unsigned char)1);

  int PosWord;
  if (SignWord<0)
    PosWord = - SignWord;
  else
    PosWord = SignWord;

  WriteBits((unsigned int)PosWord, NbBits);
}


//  Flush the buffer to the output stream
void WBitStream::Flush()
{
  pStream->flush();
}



// These are the functions specific to RBitStream

// Default constructor
RBitStream::RBitStream()
    :BitStream()
{
  piStream = 0;
}


// Constructor with the file name to open for writing
RBitStream::RBitStream(const char* pFileName)
 :BitStream()
{
  // Open the stream for writing
  piStream = new ifstream(pFileName,ios::binary);
  // Initializes position
  piStream->seekg(0);
  vFillBuffer();
}


// Destructor
RBitStream::~RBitStream()
{
  delete piStream;
}


void RBitStream::LoadValue(unsigned int NbBitsUsed, unsigned int & ValueToRead)
{
  int b_Error = Eof();
  if (b_Error)
    cout << "Error in loading" << endl;
  else
    ReadBits(ValueToRead,NbBitsUsed);
}


// Detect the end of the file
int RBitStream::Eof()
{
  return piStream->eof();
}


// Reading a bit and return it
unsigned int RBitStream::ReadBit()
{
  unsigned int Bit = ReadCurrentBit();
  if (ToNextBit()) vFillBuffer();
  return Bit;
}


// reading a bit and check the end of file
int RBitStream::ReadBit(unsigned int& Bit)
{
  Bit = ReadCurrentBit();
  if (ToNextBit())
  {
    if (!vFillBuffer()) return (0);
  }
  return (1);
}


// Read bits from the file
// Warning you must have : ui_NbBits <= sizeof(unsigned int)*8
unsigned int RBitStream::ReadBits(unsigned int NbBits)
{
  unsigned int Word;
  unsigned int EffectNbBits;

  EffectNbBits= ReadBits(Word,NbBits);

  if (EffectNbBits != NbBits)
  {
    cout << "ReadBits: Nb bits effectively read: " << EffectNbBits << endl;
  }

  return Word;
}


// Read bits from the file
float RBitStream::ReadFloat()
{
  unsigned int ToRead;
  ReadBits(ToRead,sizeof(unsigned int)*8);
  float* pValue = (float*) &ToRead;
  return (*pValue);
}


// Read bits from the file
double RBitStream::ReadDouble()
{
  unsigned int ToRead[2];
  ReadBits(ToRead[0],sizeof(unsigned int)*8);
  ReadBits(ToRead[1],sizeof(unsigned int)*8);

  double* pValue = (double*) ToRead;
  return (*pValue);
}


// Read bits from the file
unsigned char RBitStream::ReadChar()
{
  unsigned char Value;

  if (BitPos < 7)
  {
    while (!ToNextBit());
    vFillBuffer();
    Value = Buffer;
    vFillBuffer();
  }
  else
  {
    Value = Buffer;
    vFillBuffer();
  }

  Counter += 8;
  return Value;
}


// Read bits from the file
int RBitStream::ReadSignBits(unsigned int NbBits)
{
  int Value;
  int Sign;

  if (ReadBit() == 0) Sign = 1;
  else Sign = -1;

  Value = ReadBits(NbBits);

  return Value * Sign;
}


// Read bits from the file
// Warning we must have : NbBits <= sizeof(unsigned int)*8)
unsigned int RBitStream::ReadBits(unsigned int & Word, unsigned int NbBits)
{
  unsigned int WordMask = 1;
  unsigned int BitPos;
  unsigned int BitVal;
  unsigned int GoOn = (NbBits > 0);

  Word = 0;
  BitPos = 0;
  while (GoOn) {
    GoOn = ReadBit(BitVal);
    // the reading is valid
    if (BitVal == 1) Word|=WordMask;
    WordMask <<= 1;
    BitPos++;
    if (BitPos >= NbBits) GoOn = 0;
  }

  return BitPos;
}


// Fill the buffer
unsigned int RBitStream::vFillBuffer()
{
  if (!(piStream->get(Buffer)))
  {
    Buffer = 0;
    return 0;
  }
  return 1;
}
