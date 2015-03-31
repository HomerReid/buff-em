#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <signal.h>
#include <unistd.h>
#include <unistd.h>
#include <math.h>

#include "libhrutil.h"
#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#if defined(_WIN32)
#  include <windows.h>
#  include <process.h>
#else
#  include <sys/resource.h>
#  include <sys/times.h>
#endif

#if defined(HAVE_PTHREAD)
#  include <pthread.h>
#endif 

#if defined(USE_OPENMP)
#  include <omp.h>
#endif 

#include <libhrutil.h>
#include <libscuff.h>
#include <rwlock.h>

void InitTaskTiming(const char **pTaskNames);
void ResetTaskTiming();
void AddTaskTiming(int WhichTask, double Elapsed);
void LogTaskTiming();

using namespace scuff;

namespace buff {

/***************************************************************/
/***************************************************************/
/***************************************************************/
static int NumTasks=0;
static int *TaskCounts=0;
static double *TaskTimes=0;
static double *TaskTimes2=0;
static char **TaskNames=0;
static rwlock TaskLock;
void InitTaskTiming(const char **pTaskNames)
{
  if (TaskCounts)  free(TaskCounts);
  if (TaskTimes)   free(TaskTimes);
  if (TaskTimes2)  free(TaskTimes2);
  if (TaskNames)
   { for(int nt=0; nt<NumTasks; nt++)
      free(TaskNames[nt]);
     free(TaskNames);
   };
  
  NumTasks=0;
  while ( pTaskNames[NumTasks] )
   NumTasks++;

  TaskCounts = (int *)   mallocEC(NumTasks*sizeof(int));
  TaskTimes  = (double *)mallocEC(NumTasks*sizeof(double));
  TaskTimes2 = (double *)mallocEC(NumTasks*sizeof(double));
  TaskNames  = (char **) mallocEC(NumTasks*sizeof(char *));

  for(int nt=0; nt<NumTasks; nt++)
   TaskNames[nt]=strdup(pTaskNames[nt]);

}

void ResetTaskTiming()
{
  memset(TaskCounts, 0, NumTasks*sizeof(int));
  memset(TaskTimes , 0, NumTasks*sizeof(double));
  memset(TaskTimes2, 0, NumTasks*sizeof(double));
}

void AddTaskTiming(int WhichTask, double Elapsed)
{
  if (WhichTask>=NumTasks)
   { //Warn("%s:%i: internal error (%i>%i)",__FILE__,__LINE__,WhichTask,NumTasks);  
     return;
   }
  TaskLock.write_lock();
  TaskCounts[WhichTask]++;
  TaskTimes[WhichTask]+=Elapsed;
  TaskTimes2[WhichTask]+=Elapsed*Elapsed;
  TaskLock.write_unlock();
}

void LogTaskTiming(const char *Title)
{
  if (Title) Log("Task tally for %s: ",Title);
  for(int nt=0; nt<NumTasks; nt++)
   { 
     int N = TaskCounts[nt];
     if (N==0) N=1;
     double Avg    = TaskTimes[nt]  / N;
     double Avg2   = TaskTimes2[nt] / N;
     double StdDev = sqrt(Avg2 - Avg*Avg);
     Log("Task %s: (Count, Total, Mean, StdDev)={%8i,%.2e,%.2e,%.2e}",
          TaskNames[nt], TaskCounts[nt], TaskTimes[nt], Avg, StdDev);
   };
  
}

} // namespace buff
