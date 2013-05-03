#include "protein.h"

void Protein::extract_pos(size_t pos, size_t &pack_index, size_t &bit_index) {
  // TODO: look-up table to get pack index and bit index
  pos *= AA_BITS;
  pack_index = pos >> 6; // pos / PACK_BITS = 2^6
  bit_index = pos - (pack_index << 6); // pos % PACK_BITS
}

void Protein::push(AminoAcid value) {
  if (packs_.size() == 0) {
    packs_.push_back(value);
  }

  pack_t& pack = packs_.back();
  //pack &= ~mask(bit_index_, AA_BITS);  // clear current value, not really neccessary here but just to be safe
  pack &= ~(AA_MASK << bit_index_);
  pack |= ((pack_t)value) << bit_index_;
  bit_index_ += AA_BITS;
  if (bit_index_ >= PACK_BITS) {
    bit_index_ -= PACK_BITS;
    //printf("bit_index=%d, remain=%d\n", bit_index_, value>>(AA_BITS-bit_index_));
    packs_.push_back(value >> (AA_BITS - bit_index_));
  }
  ++size_;
}

void Protein::set(size_t pos, AminoAcid value) {
  size_t pack_index, bit_index;
  extract_pos(pos, pack_index, bit_index);

  pack_t& pack = packs_[pack_index];
  pack &= ~mask(bit_index, AA_BITS);  // clear current value
  pack |= ((pack_t)value) << bit_index;  // set new value

  int saved_bits = PACK_BITS - bit_index;
  int remain_bits = AA_BITS - saved_bits;
  if (remain_bits > 0) {
    pack_t& next_pack = packs_[pack_index + 1];
    next_pack &= PACK_ONES << remain_bits;  // clear 
    next_pack |= value >> saved_bits;  // set
  }
}

AminoAcid Protein::get(size_t pos) const {
  // get the index of pack element in packs_ vector that contains this AA
  // position, as well as the bit index in that pack
  size_t pack_index, bit_index;
  extract_pos(pos, pack_index, bit_index);

  // PACK_BITS = 8
  // bit_index = 5
  // pack[k]   = 11100110
  //               ^
  // pack[k+1] = 10101011
  // result = 11100110 & 11100000 >> 5
  //                       
  //          11100000 >> 5
  //          00000111
  AminoAcid result = (packs_[pack_index] >> bit_index) & AA_MASK;

  int bits_extracted = PACK_BITS - bit_index; // 3
  int remain_bits = AA_BITS - bits_extracted; // 2

  if (remain_bits > 0) {
    // result = 00000111 | (10101011 & (11111111 >> 6)
    //                                  00000011
    //                       00000011 << 3
    //                       00011000
    //          00011111
    pack_t pack_mask = PACK_ONES >> (PACK_BITS - remain_bits);
    result |= (packs_[pack_index + 1] & pack_mask) << bits_extracted;
  }

  return result;
}

AminoAcid Protein::operator[](size_t pos) const {
  return get(pos);
}
