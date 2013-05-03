#include "protein.h"
#include <iostream>
#include <sstream>
#include <assert.h>
#include <bitset>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

using namespace std;

TEST(Common, Constants) {
  ASSERT_EQ(64, PACK_BITS);
  ASSERT_EQ(5, AA_BITS);
  ASSERT_EQ(0x1F, AA_MASK);
}

TEST(Mask, Basic) {
  EXPECT_EQ(0x1F, Protein::mask(0, 5));
  EXPECT_EQ(0x3e0, Protein::mask(5, 5));
  EXPECT_EQ(0x7C00, Protein::mask(10, 5));
}

TEST(Mask, Overflow) {
  EXPECT_EQ(0x8000000000000000, Protein::mask(63, 3));
}

TEST(ExtractBits, Basic) {
  EXPECT_EQ(0x0f, Protein::extract_bits(PACK_ONES, 5, 4));
}

TEST(ExtractBits, Overflow) {
  EXPECT_EQ(1, Protein::extract_bits(PACK_ONES, 63, 10));
}

class TestableProtein : public Protein {
  public:
    const std::vector<pack_t>& packs() const { return Protein::packs(); }
    const size_t bitIndex() const { return Protein::bitIndex(); }
};

void reverse_string(string& s) {
  for (int i = 0; i < s.length()/2; i++) {
    char c = s[i];
    s[i] = s[s.length()-1-i];
    s[s.length()-1-i] = c;
  }
}

//void printProteinBits(const TestableProtein& p, ostream* os) {
ostream& operator<< (ostream& os, const TestableProtein& p) {
  stringstream ss;
  for (int i=p.packs().size()-1; i>=0; i--)
    //ss << std::hex << p.packs()[i] << std::dec;
    ss << bitset<PACK_BITS>(p.packs()[i]);
  string s = ss.str();
  reverse_string(s);
  // add spaces
  string s1 = "";
  for (int i=0; i<s.length(); i++) {
    s1.push_back(s[i]);
    if (i%PACK_BITS == PACK_BITS-1)
      s1.push_back(',');
    else if (i%5==4)
      s1.push_back(' ');
  }
  reverse_string(s1);

  os << "0x" << s1;
  return os;
}

TEST(Protein, Size0WhenNew) {
  Protein p;
  EXPECT_EQ(0, p.size());
}

TEST(Protein, SizeIncreasesWhenPush) {
  Protein p;
  p.push(0);
  EXPECT_EQ(1, p.size());
  p.push(0);
  EXPECT_EQ(2, p.size());
  p.push(0);
  p.push(0);
  p.push(0);
  EXPECT_EQ(5, p.size());
}

TEST(Protein_Big, SizeIncreasesWhenPush) {
  Protein p;
  for (int i=0; i<65536; i++) {
    p.push(0);
    ASSERT_EQ(i+1, p.size());
  }
}

TEST(Protein, BitIndexInCreasesWhenPush) {
  TestableProtein p;
  p.push(1); EXPECT_EQ(5, p.bitIndex());
  p.push(1); EXPECT_EQ(10, p.bitIndex());
  p.push(1); EXPECT_EQ(15, p.bitIndex());
  p.push(1); EXPECT_EQ(20, p.bitIndex());
  p.push(1); EXPECT_EQ(25, p.bitIndex());
  p.push(1); EXPECT_EQ(30, p.bitIndex());
  p.push(1); EXPECT_EQ(35, p.bitIndex());
  p.push(1); EXPECT_EQ(40, p.bitIndex());
  p.push(1); EXPECT_EQ(45, p.bitIndex());
  p.push(1); EXPECT_EQ(50, p.bitIndex());
  p.push(1); EXPECT_EQ(55, p.bitIndex());
  p.push(1); EXPECT_EQ(60, p.bitIndex());
  p.push(1); EXPECT_EQ(1, p.bitIndex());
  p.push(1); EXPECT_EQ(6, p.bitIndex());
}

TEST(Protein_Big, BitIndexIncreaseWhenPush) {
  TestableProtein p;
  for (int i=0; i<65536; i++) {
    p.push(i % 20);
    ASSERT_EQ((i+1)*5%PACK_BITS, p.bitIndex());
  }
}

TEST(Protein, Push) {
  TestableProtein p;
  p.push(1); ASSERT_EQ(0x0000000000000001, p.packs().at(0));  //                   0 0001
  p.push(1); ASSERT_EQ(0x0000000000000021, p.packs().at(0));  //             00 0010 0001
  p.push(1); ASSERT_EQ(0x0000000000000421, p.packs().at(0));  //       000 0100 0010 0001
  p.push(1); ASSERT_EQ(0x0000000000008421, p.packs().at(0));  // 0000 1000 0100 0010 0001
  p.push(1); ASSERT_EQ(0x0000000000108421, p.packs().at(0));  
  p.push(1); ASSERT_EQ(0x0000000002108421, p.packs().at(0));  
  p.push(1); ASSERT_EQ(0x0000000042108421, p.packs().at(0));  
  p.push(1); ASSERT_EQ(0x0000000842108421, p.packs().at(0));  
  p.push(1); ASSERT_EQ(0x0000010842108421, p.packs().at(0));  
  p.push(1); ASSERT_EQ(0x0000210842108421, p.packs().at(0));  
}

TEST(Protein, Push1) {
  TestableProtein p;
  p.push(0);  ASSERT_EQ(0x0000000000000000, p.packs().at(0)); // 0 0000
  p.push(1);  ASSERT_EQ(0x0000000000000020, p.packs().at(0)); // 00 0010
  p.push(2);  ASSERT_EQ(0x0000000000000820, p.packs().at(0)); // 000 1000
  p.push(3);  ASSERT_EQ(0x0000000000018820, p.packs().at(0)); // 0001 1000
  p.push(4);  ASSERT_EQ(0x0000000000418820, p.packs().at(0)); // 0 0100
  p.push(5);  ASSERT_EQ(0x000000000A418820, p.packs().at(0)); // 00 1010
  p.push(6);  ASSERT_EQ(0x000000018A418820, p.packs().at(0)); // 001 1000
  p.push(7);  ASSERT_EQ(0x000000398A418820, p.packs().at(0)); // 0011 1001
  p.push(8);  ASSERT_EQ(0x000008398A418820, p.packs().at(0)); // 0 1000
  p.push(9);  ASSERT_EQ(0x000128398A418820, p.packs().at(0)); // 01 0010
  p.push(10); ASSERT_EQ(0x002928398A418820, p.packs().at(0)); // 010 1001
  p.push(11); ASSERT_EQ(0x05A928398A418820, p.packs().at(0)); // 0101 1010
  p.push(12); ASSERT_EQ(0xC5A928398A418820, p.packs().at(0)); // 0 1100
              ASSERT_EQ(0x0000000000000000, p.packs().at(1));
  p.push(13); ASSERT_EQ(0xC5A928398A418820, p.packs().at(0));
              ASSERT_EQ(0x000000000000001A, p.packs().at(1)); // 01 1010
  p.push(14); ASSERT_EQ(0xC5A928398A418820, p.packs().at(0));
              ASSERT_EQ(0x000000000000039A, p.packs().at(1)); // 011 1001
  p.push(15); ASSERT_EQ(0xC5A928398A418820, p.packs().at(0));
              ASSERT_EQ(0x0000000000007B9A, p.packs().at(1)); // 0111 1011
  p.push(16); ASSERT_EQ(0xC5A928398A418820, p.packs().at(0));
              ASSERT_EQ(0x0000000000107B9A, p.packs().at(1)); // 1 0000
}

/*
TEST(Protein_Big, Push) {
  TestableProtein p;
  int bit_index = 0;
  pack_t pack = 0;
  for (int i=0; i<65536; i++) {
    p.push(i%20);
    //pack += ((pack_t)(i%20)) * pow(2.0, bit_index);
    pack |= (pack_t)(i%20) << bit_index;
    cout << "i="<<i << " 0x" << hex << pack << dec << endl;
    bit_index += AA_BITS;
    if (bit_index >= PACK_BITS) {
      bit_index -= PACK_BITS;
      pack = 0;
    }
    ASSERT_EQ(pack, p.packs().back())
      << " i=" << i << " bit_index=" << bit_index
      << "\nexpected_pack=\n0x"<<hex<<pack
      << "\npack=\n0x"<<p.packs().back()
      << dec<<endl;
  }
}
*/

TEST(Protein, PushUntilFullPack) {
  EXPECT_EQ(12, PACK_BITS/AA_BITS);

  TestableProtein p;
  for (int i=0; i<PACK_BITS/AA_BITS; i++) {  // 64/5=12 pushes
    p.push(1);
    EXPECT_EQ(1, p.packs().size());
  }
  for (int i=0; i<(PACK_BITS + 4)/AA_BITS; i++) {  // still have 4 bits from prev pack
    p.push(2);
    EXPECT_EQ(2, p.packs().size());
  }
  p.push(3);
  EXPECT_EQ(3, p.packs().size());
}

#define ASSERT_POS(pos, exp_pack_index, exp_bit_index) { \
  size_t pack_index, bit_index; \
  Protein::extract_pos(pos, pack_index, bit_index); \
  ASSERT_EQ(exp_pack_index, pack_index) << "for pos="<<pos; \
  ASSERT_EQ(exp_bit_index, bit_index); \
}

TEST(Protein, ExtractPos) {
  ASSERT_POS(0, 0, 0);
  ASSERT_POS(1, 0, 5);
  ASSERT_POS(2, 0, 10);
  ASSERT_POS(3, 0, 15);
  ASSERT_POS(4, 0, 20);
  ASSERT_POS(5, 0, 25);
  ASSERT_POS(6, 0, 30);
  ASSERT_POS(7, 0, 35);
  ASSERT_POS(8, 0, 40);
  ASSERT_POS(9, 0, 45);
  ASSERT_POS(10, 0, 50);
  ASSERT_POS(11, 0, 55);
  ASSERT_POS(12, 0, 60);
  ASSERT_POS(13, 1, 1);
  ASSERT_POS(14, 1, 6);
  ASSERT_POS(15, 1, 11);
  ASSERT_POS(16, 1, 16);
}

TEST(Protein_Big, ExtractPos) {
  int exp_pack = 0, exp_bit = 0;
  for (int i=0; i<65536; i++) {
    exp_pack = i*AA_BITS / PACK_BITS;
    exp_bit = i*AA_BITS % PACK_BITS;
    ASSERT_POS(i, exp_pack, exp_bit);
  }
}

TEST(Protein, Set) {
  TestableProtein p;
  p.push(1);
  EXPECT_EQ(0x0000000000000001, p.packs().at(0));
  p.set(0, 15);
  EXPECT_EQ(0x000000000000000F, p.packs().at(0));
  p.push(5);
  EXPECT_EQ(0x00000000000000AF, p.packs().at(0));
  p.set(1, 15);
  EXPECT_EQ(0x00000000000001EF, p.packs().at(0));
}

TEST(Protein, PushThenGet) {
  Protein p;
  p.push(0); ASSERT_EQ(0, p.get(0));
  p.push(1); ASSERT_EQ(1, p.get(1));
  p.push(2); ASSERT_EQ(2, p.get(2));
  p.push(3); ASSERT_EQ(3, p.get(3));
  p.push(4); ASSERT_EQ(4, p.get(4));
  p.push(5); ASSERT_EQ(5, p.get(5));
  p.push(6); ASSERT_EQ(6, p.get(6));
  p.push(7); ASSERT_EQ(7, p.get(7));
  p.push(8); ASSERT_EQ(8, p.get(8));
  p.push(9); ASSERT_EQ(9, p.get(9));
  p.push(10); ASSERT_EQ(10, p.get(10));
  p.push(11); ASSERT_EQ(11, p.get(11));
  p.push(12); ASSERT_EQ(12, p.get(12));
  p.push(13); ASSERT_EQ(13, p.get(13));
  p.push(14); ASSERT_EQ(14, p.get(14));
  p.push(15); ASSERT_EQ(15, p.get(15));
  p.push(16); ASSERT_EQ(16, p.get(16));
}

TEST(Protein_Big, PushThenGet) {
  TestableProtein p;
  for (int i=0;i<65536; i++) {
    p.push(i%20);
    //cout << "i = " << i << ", " << p << endl;
    ASSERT_EQ(i%20, p.get(i)) << "i = " << i << " packs = " << p;
  }
  for (int i=0; i<65536; i++) {
    ASSERT_EQ(i%20, p.get(i)) << "i = " << i;
  }
}

TEST(Protein_Big, PushThenSetThenGet) {
  TestableProtein p;
  srand(time(NULL));
  for (int i=0; i<65536; i++)
    p.push(rand()*20);
  for (int i=0; i<65536; i++)
    p.set(i, i%20);
  for (int i=0; i<65536; i++)
    ASSERT_EQ(i%20, p.get(i)) << "i = " << i;
}

int main(int argc, char **argv) {
  printf("sizeof(pack_t) = %lu\n", sizeof(pack_t));
  printf("sizeof(AminoAcid) = %lu\n", sizeof(AminoAcid));

  ::testing::InitGoogleTest(&argc, argv);

  return RUN_ALL_TESTS();
}
