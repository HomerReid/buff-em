/* Copyright (C) 2005-2011 M. T. Homer Reid
 *
 * This file is part of BUFF-EM.
 *
 * BUFF-EM is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * BUFF-EM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

/*
 * FIBBICache.h  
 *               
 *
 * homer reid    -- 5/2015
 */

#ifndef FIBBICACHE_H
#define FIBBICACHE_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <libhrutil.h>
#include "rwlock.h"

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

using namespace scuff;

// length of a FIBBI data record in units of sizeof(double)
#define FIBBIDATALEN 6

namespace buff { 

class SWGVolume; // forward reference needed below

// routine that computes a FIBBI data record for a given basis-function pair
void ComputeFIBBIData(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      double FIBBIData[FIBBIDATALEN]);

/*--------------------------------------------------------------*/
/* 'FIBBICache' is a class that implements efficient storage    */
/* and retrieval of FIBBIData structures for many pairs of      */
/* tetrahedra.                                                  */
/*--------------------------------------------------------------*/
class FIBBICache
 { 
  public:

    // constructor, destructor 
    FIBBICache(char *MeshFileName=0);
    ~FIBBICache(); 

    // look up an entry 
    void GetFIBBIData(SWGVolume *VA, int nfA,
                      SWGVolume *VB, int nfB,
                      double Data[FIBBIDATALEN]);

    // get the number of records
    int Size();

    // store/retrieve cache to/from binary file
    void Store(const char *FileName);
    int PreLoad(const char *FileName);

    int Hits, Misses;

  private:

    // any implementation of this class will have some kind of
    // storage table, but to allow maximal flexibility in implementation
    // i am just going to store an opaque pointer to this table
    // in the class body, with all the details left up to the
    // implementation
    void *opTable;

    pthread_rwlock_t lock;

    char *PreloadFileName;
    unsigned int RecordsPreloaded;

 };

} // namespace buff

#endif // #ifndef FIBBICACHE_H
