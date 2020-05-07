#pragma once
#ifndef BLOOMFILTER_H
#define BLOOMFILTER_H
#include <bitset>
#include "mySturct.h"
using namespace std;
const int size_of_bloomset = 1024 * 1024 * 8;

class bloomFilter
{
public:
	bitset<k * size_of_codealph> zeroUnit;
	bitset<size_of_bloomset> memUnit;
	bool isInBloomFilter(string value)
	{

	}
private:
	int hashzero(string value)
	{

	}

};



#endif // !BLOOMFILTER_H

