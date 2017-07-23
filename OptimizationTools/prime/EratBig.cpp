///
/// @file   EratBig.cpp
/// @brief  Segmented sieve of Eratosthenes optimized for big sieving
///         primes. This is an optimized implementation of Tomas
///         Oliveira e Silva's cache-friendly bucket sieve algorithm:
///         http://www.ieeta.pt/~tos/software/prime_sieve.html
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "primesieve/config.hpp"
#include "primesieve/EratBig.hpp"
#include "primesieve/pmath.hpp"
#include "primesieve/primesieve_error.hpp"
#include "primesieve/Wheel.hpp"

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

using namespace std;

namespace primesieve {

/// @stop:       Upper bound for sieving
/// @sieveSize:  Sieve size in bytes
/// @limit:      Sieving primes must be <= limit,
///              usually limit = sqrt(stop)
///
EratBig::EratBig(uint64_t stop, uint_t sieveSize, uint_t limit) :
  Modulo210Wheel_t(stop, sieveSize),
  limit_(limit),
  log2SieveSize_(ilog2(sieveSize)),
  moduloSieveSize_(sieveSize - 1),
  stock_(nullptr)
{
  // '>> log2SieveSize' requires power of 2 sieveSize
  if (!isPowerOf2(sieveSize))
    throw primesieve_error("EratBig: sieveSize must be a power of 2");
  init(sieveSize);
}

void EratBig::init(uint_t sieveSize)
{
  uint_t maxSievingPrime  = limit_ / NUMBERS_PER_BYTE;
  uint_t maxNextMultiple  = maxSievingPrime * getMaxFactor() + getMaxFactor();
  uint_t maxMultipleIndex = sieveSize - 1 + maxNextMultiple;
  uint_t maxSegmentCount  = maxMultipleIndex >> log2SieveSize_;
  uint_t size = maxSegmentCount + 1;

  // EratBig uses up to 1.6 gigabytes of memory
  memory_.reserve(((1u << 30) * 2) / config::BYTES_PER_ALLOC);

  lists_.resize(size, nullptr);
  for (uint_t i = 0; i < size; i++)
    pushBucket(i);
}

/// Add a new sieving prime to EratBig
void EratBig::storeSievingPrime(uint_t prime, uint_t multipleIndex, uint_t wheelIndex)
{
  assert(prime <= limit_);
  uint_t sievingPrime = prime / NUMBERS_PER_BYTE;
  uint_t segment = multipleIndex >> log2SieveSize_;
  multipleIndex &= moduloSieveSize_;
  if (!lists_[segment]->store(sievingPrime, multipleIndex, wheelIndex))
    pushBucket(segment);
}

/// Add an empty bucket to the front of lists_[segment]
void EratBig::pushBucket(uint_t segment)
{
  // allocate new buckets
  if (!stock_)
  {
    int N = config::BYTES_PER_ALLOC / sizeof(Bucket);
    memory_.emplace_back(unique_ptr<Bucket[]>(new Bucket[N]));
    Bucket* bucket = memory_.back().get();

    for (int i = 0; i < N - 1; i++)
      bucket[i].setNext(&bucket[i + 1]);
    bucket[N-1].setNext(nullptr);
    stock_ = bucket;
  }
  Bucket* emptyBucket = stock_;
  stock_ = stock_->next();
  moveBucket(*emptyBucket, lists_[segment]);
}

void EratBig::moveBucket(Bucket& src, Bucket*& dest)
{
  src.setNext(dest);
  dest = &src;
}

/// Cross-off the multiples of big sieving
/// primes from the sieve array
///
void EratBig::crossOff(byte_t* sieve)
{
  while (lists_[0]->hasNext() || !lists_[0]->empty())
  {
    Bucket* bucket = lists_[0];
    lists_[0] = nullptr;
    pushBucket(0);
    do {
      crossOff(sieve, bucket->begin(), bucket->end());
      Bucket* processed = bucket;
      bucket = bucket->next();
      processed->reset();
      moveBucket(*processed, stock_);
    } while (bucket);
  }

  rotate(lists_.begin(), lists_.begin() + 1, lists_.end());
}

/// Cross-off the next multiple of each sieving prime in the
/// current bucket. This is an implementation of the segmented
/// sieve of Eratosthenes with wheel factorization optimized for
/// big sieving primes that have very few multiples per segment
///
void EratBig::crossOff(byte_t* sieve, SievingPrime* sPrime, SievingPrime* sEnd)
{
  Bucket** lists = &lists_[0];
  uint_t moduloSieveSize = moduloSieveSize_;
  uint_t log2SieveSize = log2SieveSize_;

  // 2 sieving primes are processed per loop iteration
  // to increase instruction level parallelism
  for (; sPrime + 2 <= sEnd; sPrime += 2)
  { 
    uint_t multipleIndex0 = sPrime[0].getMultipleIndex();
    uint_t wheelIndex0    = sPrime[0].getWheelIndex();
    uint_t sievingPrime0  = sPrime[0].getSievingPrime();
    uint_t multipleIndex1 = sPrime[1].getMultipleIndex();
    uint_t wheelIndex1    = sPrime[1].getWheelIndex();
    uint_t sievingPrime1  = sPrime[1].getSievingPrime();

    // cross-off the current multiple (unset bit)
    // and calculate the next multiple
    unsetBit(sieve, sievingPrime0, &multipleIndex0, &wheelIndex0);
    unsetBit(sieve, sievingPrime1, &multipleIndex1, &wheelIndex1);

    uint_t segment0 = multipleIndex0 >> log2SieveSize;
    uint_t segment1 = multipleIndex1 >> log2SieveSize;
    multipleIndex0 &= moduloSieveSize;
    multipleIndex1 &= moduloSieveSize;

    // move the 2 sieving primes to the list related
    // to their next multiple
    if (!lists[segment0]->store(sievingPrime0, multipleIndex0, wheelIndex0))
      pushBucket(segment0);
    if (!lists[segment1]->store(sievingPrime1, multipleIndex1, wheelIndex1))
      pushBucket(segment1);
  }

  if (sPrime != sEnd)
  {
    uint_t multipleIndex = sPrime->getMultipleIndex();
    uint_t wheelIndex    = sPrime->getWheelIndex();
    uint_t sievingPrime  = sPrime->getSievingPrime();

    unsetBit(sieve, sievingPrime, &multipleIndex, &wheelIndex);
    uint_t segment = multipleIndex >> log2SieveSize;
    multipleIndex &= moduloSieveSize;

    if (!lists[segment]->store(sievingPrime, multipleIndex, wheelIndex))
      pushBucket(segment);
  }
}

} // namespace
