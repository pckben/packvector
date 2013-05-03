#ifndef _PROTEIN_H_
#define _PROTEIN_H_

#include <vector>
#include <math.h>

typedef unsigned long long pack_t;
#define PACK_ONES (~((pack_t)0))
const int PACK_BITS = sizeof(pack_t) * 8;  // number of bits of a pack
const int AA_BITS = 5;  // number of bits of an element
const pack_t AA_MASK = PACK_ONES >> (PACK_BITS - AA_BITS);  // mask of an element
typedef char AminoAcid;

/*
enum {
  Alanine = 0, Arginine, Asparagine, AsparticAcid, Cysteine, GlutamicAcid,
  Glutamine, Glycine, Histidine, Isoleucine, Leucine, Lysine, Methionine,
  Phenylalanine, Proline, Serine, Threonine, Tryptophan, Tyrosine, Valine
};
*/

/*
*/
enum {
  Ala = 0, Arg, Asn, Asp, Cys, Glu, Gln, Gly, His, Ile, Leu, Lys, Met, Phe, Pro,
  Ser, Thr, Trp, Tyr, Val
};

/*
enum {
  A = 0, R, N, D, C, E, Q, G, H, I, L, K, M, F, P, S, T, W, Y, V
};
*/

class Protein {
 public:
   Protein() : size_(0), bit_index_(0) {}

   size_t size() const { return size_; }
   void push(AminoAcid value);
   void set(size_t pos, AminoAcid value);
   AminoAcid get(size_t pos) const;
   AminoAcid operator[] (size_t pos) const;

 public:
  static void extract_pos(size_t pos, size_t &pack_index, size_t &bit_index);
  static inline pack_t mask(size_t from_bit, size_t nbits) {
    return (PACK_ONES << from_bit) - (PACK_ONES << from_bit << nbits);
  }
  static inline pack_t extract_bits(pack_t pack, size_t from_bit, size_t nbits) {
    return (pack & mask(from_bit, nbits)) >> from_bit;
  }

 protected:
  const std::vector<pack_t>& packs() const { return packs_; }
  const size_t bitIndex() const { return bit_index_; }

 private:
   std::vector<pack_t> packs_;
   size_t size_;
   size_t bit_index_;
};

#endif  // _PROTEIN_H_
