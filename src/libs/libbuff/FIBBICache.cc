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
 * FIBBI.cc -- 
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

#include <libhrutil.h>

#include "libbuff.h"

namespace buff {

void ComputeFIBBIData(SWGVolume *OA, int nfA,
                      SWGVolume *OB, int nfB,
                      FIBBIData *Data);

#define KEYLEN 27 
#define KEYSIZE (KEYLEN*sizeof(float))

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
typedef struct
 { float Key[KEYLEN];
 } KeyStruct;

typedef std::pair<KeyStruct, FIBBIData *> KeyValuePair;

struct KeyHash
 {
   long operator() (const KeyStruct &K) const 
    { return HashFunction(K.Key); }
 };

typedef struct 
 { 
   bool operator()(const KeyStruct &K1, const KeyStruct &K2) const
    { 
      if ( memcmp( (const void *)K1.Key, 
                   (const void *)K2.Key, 
                   KEYSIZE
                 ) 
         ) return false;
      return true;
    };

 } KeyCmp;

typedef std::tr1::unordered_map< KeyStruct,
                                 FIBBIData *,
                                 KeyHash,
                                 KeyCmp> KeyValueMap;

/*--------------------------------------------------------------*/
/*- class constructor ------------------------------------------*/
/*--------------------------------------------------------------*/
FIBBICache::FIBBICache()
{
  KeyValueMap *KVM=new KeyValueMap;
  opTable = (void *)KVM;
  PreloadFileName=0;
  RecordsPreloaded=0;
  Hits=Misses=0;
}

/*--------------------------------------------------------------*/
/*- class destructor  ------------------------------------------*/
/*--------------------------------------------------------------*/
FIBBICache::~FIBBICache()
{
  if (PreloadFileName) 
   free(PreloadFileName);

  KeyValueMap *KVM=(KeyValueMap *)opTable;
  delete KVM;
} 

/***************************************************************/
/* the following three routines implement a technique for      */
/* assigning a search key to pairs of SWG basis functions in   */
/* a translation-independent (but not rotation-indepedent) way */
/***************************************************************/
static void inline VecSubFloat(double *V1, double *V2, float *V1mV2)
{ V1mV2[0] = ((float)V1[0]) - ((float)V2[0]);
  V1mV2[1] = ((float)V1[1]) - ((float)V2[1]);
  V1mV2[2] = ((float)V1[2]) - ((float)V2[2]);
}

void GetFIBBICacheKey(SWGVolume *OA, int nfA, 
                      SWGVolume *OB, int nfB,
                      KeyStruct *KS)
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

  float *Key = KS->Key;
  VecSubFloat(QMA, QPA, Key + 0*3);
  VecSubFloat(V1A, QPA, Key + 1*3);
  VecSubFloat(V2A, QPA, Key + 2*3);
  VecSubFloat(V3A, QPA, Key + 3*3);
  VecSubFloat(QPB, QPA, Key + 4*3);
  VecSubFloat(QMB, QPA, Key + 5*3);
  VecSubFloat(V1B, QPA, Key + 6*3);
  VecSubFloat(V2B, QPA, Key + 7*3);
  VecSubFloat(V3B, QPA, Key + 8*3);

}

/*--------------------------------------------------------------*/
/*- routine for fetching a FIBBI data record from a FIBBICache: */
/*- we look through our table to see if we have a record for    */
/*- this panel pair, we return it if we do, and otherwise we    */
/*- compute a new FIBBI data record for this panel pair and     */
/*- add it to the table.                                        */
/*--------------------------------------------------------------*/
FIBBIData *FIBBICache::GetFIBBIData(SWGVolume *OA, int nfA,
                                    SWGVolume *OB, int nfB)
{
  KeyStruct KS;
  GetFIBBICacheKey(OA, nfA, OB, nfB, &KS);

  /***************************************************************/
  /* look for this key in the cache ******************************/
  /***************************************************************/
  KeyValueMap *KVM=(KeyValueMap *)opTable;

  FCLock.read_lock();
  KeyValueMap::iterator p=KVM->find(KS);
  FCLock.read_unlock();

  if ( p != (KVM->end()) )
   { Hits++;
     return (FIBBIData *)(p->second);
   };
  
  /***************************************************************/
  /* if it was not found, allocate and compute a new FIBBIData   */
  /* structure, then add this structure to the cache             */
  /***************************************************************/
  Misses++;
  KeyStruct *K2 = (KeyStruct *)mallocEC(sizeof(*K2));
  memcpy(K2->Key, KS.Key, KEYSIZE);
  FIBBIData *Data=(FIBBIData *)mallocEC(sizeof(FIBBIData));
  ComputeFIBBIData(OA, nfA, OB, nfB, Data);
   
  FCLock.write_lock();
  KVM->insert( KeyValuePair(*K2, Data) );
  FCLock.write_unlock();

  return Data;
}

/***************************************************************/
/* this routine and the following routine implement a mechanism*/
/* for storing the contents of a FIBBI cache to a binary file, */
/* and subsequently pre-loading a FIBBI cache with the content */
/* of a file created by this storage operation.                */
/*                                                             */
/* the file format is pretty simple (and non-portable w.r.t.   */
/* endianness):                                                */
/*  bytes 0--10:   'FIBBICACHE' + 0 (a file signature used as  */
/*                                   a simple sanity check)    */
/*  next xx bytes:  first record                               */
/*  next xx bytes:  second record                              */
/*  ...             ...                                        */
/*                                                             */
/* where xx is the size of the record; each record consists of */
/* a search key (KEYLEN float values) followed by the content  */
/* of the FIBBIDataRecord for that search key.                 */
/*                                                             */
/* note: FIBBICF = 'FIBBI cache file'                          */
/***************************************************************/
const char FIBBICF_Signature[]="FIBBICACHE";
#define FIBBICF_SIGSIZE sizeof(FIBBICF_Signature)

// note that this structure differs from the KeyValuePair structure 
// defined above in that it contains the actual contents of 
// QIFIBBIData structure, whereas KeyValuePair contains just a 
// pointer to such a structure.
typedef struct FIBBICF_Record
 { KeyStruct K;
   FIBBIData DataBuffer;
 } FIBBICF_Record;
#define FIBBICF_RECSIZE sizeof(FIBBICF_Record)

void FIBBICache::Store(const char *FileName)
{
  KeyValueMap *KVM=(KeyValueMap *)opTable;

  if (FileName==0) return;

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
  /*--------------------------------------------------------------*/
  /*--------------------------------------------------------------*/
  FCLock.read_lock();

  KeyValueMap::iterator it;
  int NumRecords=0;

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(FileName,"w");
  if (!f)
   { fprintf(stderr,"warning: could not open file %s (aborting cache dump)...",FileName);
     goto done;
   };
  Log("Writing FIBBI cache to file %s...",FileName);

  /*---------------------------------------------------------------------*/
  /*- write file signature ----------------------------------------------*/
  /*---------------------------------------------------------------------*/
  fwrite(FIBBICF_Signature,FIBBICF_SIGSIZE,1,f);

  /*---------------------------------------------------------------------*/
  /*- iterate through the table and write entries to the file one-by-one */
  /*---------------------------------------------------------------------*/
  KeyStruct K;
  FIBBIData *Data;
  FIBBICF_Record MyRecord;
  for ( it = KVM->begin(); it != KVM->end(); it++ ) 
   { 
     K=it->first;
     Data=it->second;
     memcpy(&(MyRecord.K.Key),      K.Key, sizeof(MyRecord.K.Key ));
     memcpy(&(MyRecord.DataBuffer), Data,  sizeof(MyRecord.DataBuffer));
     if ( 1 != fwrite(&(MyRecord),sizeof(MyRecord),1,f ) )
      break;
     NumRecords++;
   };

  /*---------------------------------------------------------------------*/
  /*- and that's it -----------------------------------------------------*/
  /*---------------------------------------------------------------------*/
  fclose(f);
  Log(" ...wrote %i FIBBI records.",NumRecords);

 done:
  FCLock.read_unlock();
}

void FIBBICache::PreLoad(const char *FileName)
{

  FCLock.write_lock();

  KeyValueMap *KVM=(KeyValueMap *)opTable;

  /*--------------------------------------------------------------*/
  /*- try to open the file ---------------------------------------*/
  /*--------------------------------------------------------------*/
  FILE *f=fopen(FileName,"r");
  if (!f)
   { fprintf(stderr,"warning: could not open file %s (skipping cache preload)\n",FileName);
     Log("Could not open FIBBI cache file %s...",FileName); 
     goto done;
   };

  /*--------------------------------------------------------------*/
  /*- run through some sanity checks to make sure we have a valid */
  /*- cache file                                                  */
  /*--------------------------------------------------------------*/
  const char *ErrMsg;
  ErrMsg = 0;
  
  struct stat fileStats;
  if ( fstat(fileno(f), &fileStats) )
   ErrMsg="invalid cache file";
  
  // check that the file signature is present and correct 
  off_t FileSize;
  FileSize=fileStats.st_size;
  char FileSignature[FIBBICF_SIGSIZE]; 
  if ( ErrMsg==0 && (FileSize < (signed int )FIBBICF_SIGSIZE) )
   ErrMsg="invalid cache file";
  if ( ErrMsg==0 && 1!=fread(FileSignature, FIBBICF_SIGSIZE, 1, f) )
   ErrMsg="invalid cache file";
  if ( ErrMsg==0 && strcmp(FileSignature, FIBBICF_Signature) ) 
   ErrMsg="invalid cache file";

  // the file size, minus the portion taken up by the signature, 
  // should be an integer multiple of the size of a FIBBICF_Record
  FileSize-=FIBBICF_SIGSIZE;
  if ( ErrMsg==0 && (FileSize % FIBBICF_RECSIZE)!=0 )
   ErrMsg="cache file has incorrect size";

  /*--------------------------------------------------------------*/
  /*- allocate memory to store cache entries read from the file. -*/
  /*- for now, we abort if there is not enough memory to store   -*/
  /*- the full cache; TODO explore schemes for partial preload.  -*/
  /*--------------------------------------------------------------*/
  unsigned int NumRecords;
  FIBBICF_Record *Records;
  NumRecords = FileSize / FIBBICF_RECSIZE;
  if ( ErrMsg==0 )
   { Records=(FIBBICF_Record *)mallocEC(NumRecords*FIBBICF_RECSIZE);
     if ( !Records)
      ErrMsg="insufficient memory to preload cache";
   };

  /*--------------------------------------------------------------*/
  /*- pause here to clean up if anything went wrong --------------*/
  /*--------------------------------------------------------------*/
  if (ErrMsg)
   { fprintf(stderr,"warning: file %s: %s (skipping cache preload)\n",FileName,ErrMsg);
     Log("FIBBI cache file %s: %s (skipping cache preload)",FileName,ErrMsg);
     fclose(f);
     goto done;
   };

  /*--------------------------------------------------------------*/
  /*- now just read records from the file one at a time and add   */
  /*- them to the table.                                          */
  /*--------------------------------------------------------------*/
  unsigned int nr;
  Log("Preloading FIBBI records from file %s...",FileName);
  for(nr=0; nr<NumRecords; nr++)
   { 
     if ( fread(Records+nr, FIBBICF_RECSIZE,1,f) != 1 )
      { fprintf(stderr,"warning: file %s: read only %i of %i records",FileName,nr+1,NumRecords);
        fclose(f);
	goto done;
      };

     KVM->insert( KeyValuePair(Records[nr].K, &(Records[nr].DataBuffer)) );
   };

  /*--------------------------------------------------------------*/
  /*- the full file was successfully preloaded -------------------*/
  /*--------------------------------------------------------------*/
  Log(" ...successfully preloaded %i FIBBI records.",NumRecords);
  fclose(f);

  // the most recent file from which we preloaded, and the number of 
  // records preloaded, are stored within the class body to allow us 
  // to skip dumping the cache back to disk in cases where that would
  // amount to just rewriting the same cache file
  if (PreloadFileName)
   free(PreloadFileName);
  PreloadFileName=strdupEC(FileName);
  RecordsPreloaded=NumRecords;

 done:
  FCLock.write_unlock();
}

} // namespace buff
