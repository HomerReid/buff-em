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
 * FIBBICache.cc -- caching of (F)requency-(I)ndependent 
 *               -- (B)asis function -- (B)asis function (I)ntegrals
 * 
 * homer reid    -- 4/2015
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <tr1/unordered_map>
#include <pthread.h> // needed for rwlock
#include <omp.h> // needed for rwlock

#include <libhrutil.h>

#include "libbuff.h"

namespace buff {

void ComputeFIBBIData(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      double *GFI, double *dGFI);

#define KEYLEN 27
#define KEYSIZE (KEYLEN*sizeof(float))

// for GCache entries, only 6 of the 48 entries in Data are used.
typedef struct DataRecord
 { float Key[27];
   double Data[48];
 } DataRecord;

#define GCACHE_DATALEN 6
#define DGCACHE_DATALEN 48

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*- note: i found this on wikipedia ... ------------------------*/
/*--------------------------------------------------------------*/
static long JenkinsHash(const char *key, size_t len)
{
    long hash;
    unsigned int i;
    for(hash = i = 0; i < len; ++i)
    {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}

long HashFunction(const float *Key)
{ return JenkinsHash( (const char *)Key, KEYSIZE );
} 

/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
/*--------------------------------------------------------------*/
typedef float  *KeyType;
typedef double *DataType;

typedef std::pair<KeyType, DataType> KeyValuePair;

struct KeyHash
 {
   long operator() (const KeyType &K) const 
    { return HashFunction(K); }
 };

typedef struct 
 { 
   bool operator()(const KeyType &K1, const KeyType &K2) const
    { 
      return !memcmp( (const void *)K1, (const void *)K2, KEYSIZE );
    };

 } KeyCmp;

typedef std::tr1::unordered_map< KeyType,
                                 DataType,
                                 KeyHash,
                                 KeyCmp> KeyValueMap;

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIBBICache::FIBBICache(char *MeshFileName, bool pIsGCache)
{
  KeyValueMap *KVM=new KeyValueMap;
  opTable = (void *)KVM;
  IsGCache=pIsGCache;

  pthread_rwlock_init(&lock,0);

  Hits=Misses=0;
  RecordBuffer=0;
  NumRecords=RecordBufferLen=0;
  RBChunkSize=1000; // default chunk size for reallocation
  char *s=getenv("BUFF_FIBBI_CHUNK");
  //if (char *s=getenv("BUFF_FIBBI_CHUNK"))
  if (s)
   { sscanf(s,"%i",&RBChunkSize);
     Log("Set FIBBI chunk size to %i.",RBChunkSize);
   };

  PreloadFileName=0;
  RecordsPreloaded=0;
  if (MeshFileName)
   { char CacheFileName[MAXSTR];
     const char *Extension = IsGCache ? "GCache" : "dGCache";
     snprintf(CacheFileName, MAXSTR, "%s.%s", 
              RemoveExtension(MeshFileName), Extension);
     PreLoad(CacheFileName);
   };
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIBBICache::~FIBBICache()
{
  pthread_rwlock_destroy(&lock);

  if (PreloadFileName)
   free(PreloadFileName);

  if (RecordBuffer)
   free(RecordBuffer);

  KeyValueMap *KVM=(KeyValueMap *)opTable;
  delete KVM;

} 

/***************************************************************/
/* the following three routines implement a technique for      */
/* assigning a search key to pairs of SWG basis functions in   */
/* a translation-independent (but not rotation-indepedent) way */
/***************************************************************/
static void inline VecSubFloat(double *V1, double *V2, float *V1mV2)
{ V1mV2[0] = (float)(V1[0] - V2[0]);
  V1mV2[1] = (float)(V1[1] - V2[1]);
  V1mV2[2] = (float)(V1[2] - V2[2]);
}

void GetFIBBICacheKey(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      KeyType K)
{

  double  *VA = OA->Vertices;
  SWGFace *FA = OA->Faces[nfA];
  double *QPA = VA + 3*(FA->iQP);
  double *QMA = VA + 3*(FA->iQM);
  double *V1A = VA + 3*(FA->iV1);
  double *V2A = VA + 3*(FA->iV2);
  double *V3A = VA + 3*(FA->iV3);

  double  *VB = OB->Vertices;
  SWGFace *FB = OB->Faces[nfB];
  double *QPB = VB + 3*(FB->iQP);
  double *QMB = VB + 3*(FB->iQM);
  double *V1B = VB + 3*(FB->iV1);
  double *V2B = VB + 3*(FB->iV2);
  double *V3B = VB + 3*(FB->iV3);

  VecSubFloat(QMA, QPA, K + 0*3);
  VecSubFloat(V1A, QPA, K + 1*3);
  VecSubFloat(V2A, QPA, K + 2*3);
  VecSubFloat(V3A, QPA, K + 3*3);
  VecSubFloat(QPB, QPA, K + 4*3);
  VecSubFloat(QMB, QPA, K + 5*3);
  VecSubFloat(V1B, QPA, K + 6*3);
  VecSubFloat(V2B, QPA, K + 7*3);
  VecSubFloat(V3B, QPA, K + 8*3);

}

/*--------------------------------------------------------------*/
/*- routine for fetching a FIBBI data record from a FIBBICache: */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIBBI data record for this panel pair and     */
/*- add it to the table.                                        */
/*--------------------------------------------------------------*/
void FIBBICache::GetFIBBIData(SWGVolume *OA, int nfA,
                              SWGVolume *OB, int nfB,
                              double *Data)
{
  int DataLen    = (IsGCache) ? GCACHE_DATALEN : DGCACHE_DATALEN;
  int DataSize   = DataLen * sizeof(double);

  /***************************************************************/
  /* look for this key in the cache ******************************/
  /***************************************************************/
  float Key[KEYLEN], *K=Key;
  GetFIBBICacheKey(OA, nfA, OB, nfB, K);
  KeyValueMap *KVM=(KeyValueMap *)opTable;

  //FCLock.read_lock();
  pthread_rwlock_rdlock(&lock);
  KeyValueMap::iterator p=KVM->find(K);
  bool Found = (p != (KVM->end()) );
  if ( Found )
   { Hits++;
     memcpy(Data, (double *)(p->second), DataSize);
   };
  pthread_rwlock_unlock(&lock);
  if (Found) return;
  
  /***************************************************************/
  /* if it was not found, compute a new FIBBI data record...     */
  /***************************************************************/
  Misses++;
  
  double DataBuffer[48];
  if (IsGCache)
   ComputeFIBBIData(OA, nfA, OB, nfB, DataBuffer, 0);
  else
   ComputeFIBBIData(OA, nfA, OB, nfB, 0, DataBuffer);
  memcpy(Data, DataBuffer, DataSize);
   
  /***************************************************************/
  /* ...and add this data record to the cache                    */
  /***************************************************************/
  pthread_rwlock_wrlock(&lock);

  NumRecords++;
  int RecordSize = KEYSIZE + DataSize;
  if (NumRecords>RecordBufferLen)
   { 
     RecordBufferLen += RBChunkSize;
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
printf("thread %i reallocing to %8i \n", omp_get_thread_num(),RecordBufferLen);
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
     RecordBuffer = reallocEC(RecordBuffer, RecordBufferLen*RecordSize);
   };
  char *Record    = ((char *)RecordBuffer) + NumRecords*RecordSize;
  float *NewKey   = (float *)Record;
  float *fNewData = NewKey + KEYLEN;
  double *NewData = (double *)fNewData;

  memcpy(NewKey, Key, KEYSIZE);
  memcpy(NewData, DataBuffer, DataSize);
  KVM->insert( KeyValuePair(NewKey, NewData) );

  pthread_rwlock_unlock(&lock);
}

/***************************************************************/
/* this routine and the following routine implement a mechanism*/
/* for storing the contents of a FIBBI cache to a binary file, */
/* and subsequently pre-loading a FIBBI cache with the content */
/* of a file created by this storage operation.                */
/*                                                             */
/* the file format is pretty simple (and non-portable w.r.t.   */
/* endianness):                                                */
/*  bytes 0--12:   'FIBBI_GCACHE' + 0                          */
/*   or                                                        */
/*                 'FIBBIdGCACHE' + 0                          */
/*  next xx bytes:  first record                               */ 
/*  next xx bytes:  second record                              */
/*  ...             ...                                        */
/*                                                             */
/* where xx is the size of the record; each record consists of */
/* a search key (KEYLEN float values) followed by the content  */
/* of the FIBBI data record for that search key.               */
/*                                                             */
/* note: FIBBICF = 'FIBBI cache file'                          */
/***************************************************************/
const char FIBBICF_GSignature[]  = "FIBBI_GCACHE";
const char FIBBICF_dGSignature[] = "FIBBIdGCACHE";
#define FIBBICF_SIGSIZE (sizeof(FIBBICF_GSignature))

void FIBBICache::Store(const char *FileName)
{
  KeyValueMap *KVM=(KeyValueMap *)opTable;

  if (!FileName)
   return;

  /*--------------------------------------------------------------*/
  /*- pause to check if the following conditions are satisfied:  -*/
  /*-  (1) the FIBBI cache was preloaded from an input file whose-*/
  /*-      name matches the name of the output file to which we  -*/
  /*-      are being asked to dump the cache                     -*/
  /*-  (2) the number of cache records hasn't changed since we   -*/
  /*-      preloaded from the input file.                        -*/
  /*- if both conditions are satisfied, we don't bother to dump  -*/
  /*- the cache since the operation would result in a cache dump -*/
  /*- file identical to the one that already exists.             -*/
  /*--------------------------------------------------------------*/
  if (     PreloadFileName
       && !strcmp(PreloadFileName, FileName)
       && RecordsPreloaded==(KVM->size())
     )
   { Log("FIBBI cache unchanged since reading from %s (skipping cache dump)",FileName);
     return;
   };

  /*--------------------------------------------------------------*/
  /*- i assume that Preload() and Store() won't be called from    */
  /*- multithreaded code sections.                                */
  /*--------------------------------------------------------------*/
  //FCLock.read_lock();

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(FileName,"w");
  if (!f)
   { fprintf(stderr,"warning: could not open file %s (aborting cache dump)...",FileName);
     //FCLock.read_unlock();
     return;
   };
  Log("Writing FIBBI cache to file %s...",FileName);

  /*---------------------------------------------------------------------*/
  /*- write file signature ----------------------------------------------*/
  /*---------------------------------------------------------------------*/
  if (IsGCache)
   fwrite(FIBBICF_GSignature, FIBBICF_SIGSIZE,1,f);
  else
   fwrite(FIBBICF_dGSignature,FIBBICF_SIGSIZE,1,f);

  /*---------------------------------------------------------------------*/
  /*- iterate through the table and write entries to the file one-by-one */
  /*---------------------------------------------------------------------*/
  int NumWritten=0;
  int DataLen  = (IsGCache) ? GCACHE_DATALEN : DGCACHE_DATALEN;
  int DataSize = DataLen * sizeof(double);
  int RecordSize = KEYSIZE + DataSize;
  KeyValueMap::iterator it;
  for ( it = KVM->begin(); it != KVM->end(); it++ )
   { 
     float *K=it->first;
     if ( 1 != fwrite(K, RecordSize, 1, f ) ) break;
     NumWritten++;
   };

  /*---------------------------------------------------------------------*/
  /*- and that's it -----------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  fclose(f);
  Log(" ...wrote %i/%i FIBBI records.",NumWritten,NumRecords);

  //FCLock.read_unlock();
}

void FIBBICache::PreLoad(const char *FileName)
{
  int DataLen    = (IsGCache) ? GCACHE_DATALEN : DGCACHE_DATALEN;
  int DataSize   = DataLen * sizeof(double);
  int RecordSize = DataSize + KEYSIZE;
  KeyValueMap *KVM=(KeyValueMap *)opTable;

  //FCLock.write_lock();

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(FileName,"r");
  if (!f)
   { fprintf(stderr,"warning: could not open file %s (skipping cache preload)\n",FileName);
     Log("Could not open FIBBI cache file %s...",FileName); 
     //FCLock.write_unlock();
     return;
   };

  /*--------------------------------------------------------------*/
  /*- run through some sanity checks to make sure we have a valid */
  /*- cache file                                                  */
  /*--------------------------------------------------------------*/
  const char *ErrMsg=0;
  struct stat fileStats;
  off_t FileSize;
  if ( fstat(fileno(f), &fileStats) )
   { ErrMsg="invalid cache file";
     goto fail;
   };
  
  // check that the file signature is present and has the right size 
  FileSize=fileStats.st_size;
  char FileSignature[FIBBICF_SIGSIZE]; 
  if ( (FileSize < ((signed int )FIBBICF_SIGSIZE)) )
   { ErrMsg="invalid cache file";
     goto fail;
   };
  if ( 1!=fread(FileSignature, FIBBICF_SIGSIZE, 1, f) )
   { ErrMsg="invalid cache file";
     goto fail;
   };

  // check that the file signature is one that we recognize
  bool FileIsGCache;
  if (      !strcmp(FileSignature, FIBBICF_GSignature)  )
   FileIsGCache=true;
  else if ( !strcmp(FileSignature, FIBBICF_dGSignature) )
   FileIsGCache=false;
  else 
   { ErrMsg="invalid cache file";
     goto fail;
   };

  // check that the file type agrees with what we are
  if ( FileIsGCache!=IsGCache )
   { ErrMsg="cache file is wrong type";
     goto fail;
   };

  // the file size, minus the portion taken up by the signature,
  // should be an integer multiple of the size of a (Key,Data) pair
  FileSize-=FIBBICF_SIGSIZE;
  if ( (FileSize % RecordSize)!=0 )
   { ErrMsg="cache file has incorrect size";
     goto fail;
   };

  /*--------------------------------------------------------------*/
  /*- allocate memory to store cache entries read from the file. -*/
  /*- for now, we abort if there is not enough memory to store   -*/
  /*- the full cache; TODO explore schemes for partial preload.  -*/
  /*--------------------------------------------------------------*/
  NumRecords      = FileSize / RecordSize;
  RecordBufferLen = NumRecords;
  RecordBuffer    = (void *)mallocEC(NumRecords*DataSize);

  /*--------------------------------------------------------------*/
  /*- now just read records from the file one at a time and add   */
  /*- them to the table.                                          */
  /*--------------------------------------------------------------*/
  Log("Preloading FIBBI records from file %s...",FileName);
  for(int nr=0; nr<NumRecords; nr++)
   { 
     char *Record    = ((char *)RecordBuffer) + nr*RecordSize;
     if ( fread( Record, RecordSize, 1, f) != 1 )
      { fprintf(stderr,"warning: file %s: read only %i of %i records",
                FileName,nr+1,NumRecords);
        NumRecords=nr;        
        fclose(f);
        //FCLock.write_unlock();
	return;
      };
  
     float *NewKey   = (float *)Record;
     float *fNewData = NewKey + KEYLEN;
     double *NewData = (double *)fNewData;
     KVM->insert( KeyValuePair(NewKey, NewData) );
   };

  /*--------------------------------------------------------------*/
  /* the most recent file from which we preloaded, and the number */
  /* of records preloaded, are stored within the class body to    */
  /* allow us to skip dumping the cache back to disk in cases     */
  /* where that would amount to just rewriting the same cache     */
  /*--------------------------------------------------------------*/
  if (PreloadFileName)
   free(PreloadFileName);
  PreloadFileName=strdupEC(FileName);
  RecordsPreloaded=NumRecords;

  Log(" ...successfully preloaded %i FIBBI records.",NumRecords);
  fclose(f);
  //FCLock.write_unlock();
  return;

  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
 fail:
  fprintf(stderr,"warning: file %s: %s (skipping cache preload)\n",FileName,ErrMsg);
  Log("FIBBI cache file %s: %s (skipping cache preload)",FileName,ErrMsg);
  fclose(f);
  //FCLock.write_unlock();
  return;

}

} // namespace buff
