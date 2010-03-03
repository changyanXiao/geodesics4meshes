#ifndef _BITSTR_H_
#define _BITSTR_H_

#include <iostream>
#include <fstream>

using std::istream;
using std::ofstream;

//#define bool int
#define true 1
#define false 0


class BitStream
{
 public:

  // default constructor
  BitStream();

  // constructor with bitstream file
  BitStream(const char* pc_FileName, char mode);

  // Destructor
  virtual ~BitStream();

  // Reading current bit
  unsigned int ReadCurrentBit();

  // Writing current bit
  void WriteCurrentBit(unsigned int Bit);

  // Writing current bit (second signature)
  void WriteCurrentBit(unsigned char cBit);


  // Moving to the next buffer bit
  // and go back to zero if the end of buffer is reached
  // returns 1 if the end of buffer has been reached
  unsigned int ToNextBit();
  unsigned int ToPreviousBit();

  unsigned int GetCounter();


 protected:

  // Go back to start bit
  void tToStartBit();

  // Go back to start bit
  void tToEndBit();

  // Filling the buffer, and returns 0 on failure, 1 on success
  unsigned int vFillBuffer();

  // Data members

  // Buffer storing bit, before reading or writing
  char Buffer;

  // Current bit position in the buffer
  unsigned int BitPos;

  // Mask on the current bit
  int Mask;

  // Counter: total number of bits read or written
  unsigned int Counter;

};


class WBitStream:public BitStream
{
 public:

  // Default constructor
  WBitStream();

  // constructor with bitstream file
  WBitStream(const char* cFileName);

  // Destructor
  virtual  ~WBitStream();


  // Save value to bitstream with NbBitsUsed bits (the size is checked)
  // Return FALSE if an error occured
  int SaveValue(unsigned int NbBitsUsed, unsigned int ValueToEncode, char* pc_ErrorMessage = 0);

  // Stuff the end of buffer with a given value
  void Stuff(unsigned int ui_Bit);


  // Write a bit
   void WriteBit(unsigned int Bit);

   void WriteBit(unsigned char Bit);

  // Writing bits from an unsigned int
  void WriteBits(unsigned int ui_Word, unsigned int ui_NbBits);

  // Writing float
  void WriteFloat(float f_Value);

  // Writing double
   void WriteDouble(double d_Value);

  // Shift on next char and quick write byte
  void WriteByte(unsigned char Byte);
  void WriteSignBits(int i_SignWord, unsigned int ui_NbBits);

  void Flush();

 private:

  ofstream *pStream;
};


class RBitStream:public BitStream
{
 public:

  // Default constructor
  RBitStream();

  // constructor with bitstream file
  RBitStream(const char* cFileName);

  // Destructor
  virtual  ~RBitStream();

  int Eof();

  // Load value from bitstream with ui_NbBitdsUsed bits
  // Return FALSE if an error occured
  void LoadValue(unsigned int NbBitsUsed, unsigned int& ValueToRead);

  // Reading a bit and return it
   unsigned int ReadBit();


  // reading a bit and return 0 if the end of file is reached, 1 if
  // one can continue reading
   int ReadBit(unsigned int& Bit);


  // Reading bits and return them as an unsigned int
   unsigned int  ReadBits(unsigned int NbBits);

  // Reading bits and store them in an unsigned Int, returns the effective
  // number of read bits
   unsigned int  ReadBits(unsigned int & Word, unsigned int NbBits);

  // Reading float
   float ReadFloat();

  // Reading float
   double ReadDouble();

  // Reading 8 bits (i.e. a byte)
   unsigned char ReadChar();

  // Reading bits and return them as an int
  //   use 1 bit to encode sign
  //   use ui_NbBits to encode absolute value
  int ReadSignBits(unsigned int NbBits);

 private:

  unsigned int vFillBuffer();

  istream *piStream;
};

#endif
