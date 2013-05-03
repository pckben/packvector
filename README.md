A library to pack non-standard size numbers into a vector of normal numbers.

Eg. Store an array of 5-bit integers into a vector of 64-bit integers. The idea
is to compactly represent 20 different types of AminoAcids as 5-bit integers for
sending to GPU clusters (memory limited) to process
