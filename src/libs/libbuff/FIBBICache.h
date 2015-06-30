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

namespace buff { 

class SWGVolume; // forward reference needed below

/*--------------------------------------------------------------*/
/* 'FIBBICache' is a class that implements efficient storage    */
/* and retrieval of FIBBIData structures for many pairs of      */
/* tetrahedra.                                                  */
/*--------------------------------------------------------------*/
class FIBBICache
 { 
  public:

    // constructor, destructor 
    FIBBICache(char *MeshFileName=0, bool IsGCache=true);
    ~FIBBICache(); 

    // look up an entry 
    void GetFIBBIData(SWGVolume *VA, int nfA, SWGVolume *VB, int nfB,
                      double *Data);

    // get the number of records
    int Size();

    // store/retrieve cache to/from binary file
    void Store(const char *FileName);
    void PreLoad(const char *FileName);

    int Hits, Misses;
    bool IsGCache;

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
