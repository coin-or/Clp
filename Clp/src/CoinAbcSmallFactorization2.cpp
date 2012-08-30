/* $Id: CoinAbcSmallFactorization2.cpp 1448 2011-06-19 15:34:41Z stefan $ */
// Copyright (C) 2002, International Business Machines
// Corporation and others, Copyright (C) 2012, FasterCoin.  All Rights Reserved.
// This code is licensed under the terms of the Eclipse Public License (EPL).

#if CLP_HAS_ABC
#include "CoinAbcCommonFactorization.hpp"
#ifndef ABC_JUST_ONE_FACTORIZATION
#define CoinAbcTypeFactorization CoinAbcSmallFactorization
#define ABC_SMALL 4
#include "CoinAbcBaseFactorization.hpp"
#include "CoinAbcBaseFactorization2.cpp"
#endif
#endif
