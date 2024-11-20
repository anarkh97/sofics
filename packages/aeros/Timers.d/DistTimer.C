#include <cstdio>
#include <algorithm>
#include <Timers.d/DistTimer.h>
#include <Comm.d/Communicator.h>

DistTimer::DistTimer(int numTimers)
{
 numSubTimers = numTimers;

 timeData = new TimeData[numSubTimers];
}

TimeStats
DistTimer::getStats()
{
 return TimeStats(getMin(),getMax(),getAvg(),getTot());
}

void
DistTimer::print()
{
 TimeData tmin = getMin();
 TimeData tmax = getMax();
 TimeData tavg = getAvg();
 TimeData ttot = getTot();
 int i;
 for(i=0; i<numSubTimers; ++i)
   fprintf(stderr,"sub timer %d Time/Mem: %14.5f %14.5f\n",i+1,timeData[i].time/1000.0, (double) timeData[i].memory/(1024*1024));
 fprintf(stderr,"Timers %14.5f %14.5f %14.5f %14.5f\n",
                   tmin.time/1000.0,tavg.time/1000.0,
                   tmax.time/1000.0,ttot.time/1000.0);

 fprintf(stderr,"Memory: %14.5f %14.5f %14.5f %14.5f\n", (double) tmin.memory/(1024*1024), (double) tavg.memory/(1024*1024), (double) tmax.memory/(1024*1024), (double) ttot.memory/(1024*1024));
}

void
DistTimer::printOverAll()
{
 TimeData ttot = getTot();
 fprintf(stderr,"over all time/mem:  %14.5f %ld\n",ttot.time/1000.0, ttot.memory/(1024*1024));
}

const TimeData
DistTimer::getMin()
{
 // get local min
 long minMem = timeData[0].memory;
 double minTime   = timeData[0].time;
 int i;
 for(i=1; i<numSubTimers; ++i) {
   minMem  = std::min( minMem,  timeData[i].memory );
   minTime = std::min( minTime, timeData[i].time );
 }
   
#ifdef DISTRIBUTED
 double minMemory = (double) minMem;
 if(structCom) minMemory = structCom->globalMin(minMemory);
 minMem = (long) minMemory;
 if(structCom) minTime = structCom->globalMin(minTime);
#endif 

 TimeData ret;
 ret.memory = minMem;
 ret.time   = minTime;
 return ret;
}


TimeData
DistTimer::getMax()
{

 // get local max
 long maxMem = timeData[0].memory;
 double maxTime   = timeData[0].time;
 int i; 
 for(i=1; i<numSubTimers; ++i) {
   maxMem  = std::max(maxMem,  timeData[i].memory);
   maxTime = std::max(maxTime, timeData[i].time);
 }
 
 // communicate to get global max

#ifdef DISTRIBUTED
 double maxMemory = (double) maxMem;
 if(structCom) maxMemory  = structCom->globalMax(maxMemory);
 maxMem = (long) maxMemory;
 if(structCom) maxTime = structCom->globalMax(maxTime);
#endif 

 TimeData ret;
 ret.memory = maxMem;
 ret.time   = maxTime;
 return ret;

}

TimeData
DistTimer::getAvg()
{
 int nSub = numSubTimers;

 // get local tot
 long totMem = 0;
 double totTime = 0.0;
 int i;
 for(i=0; i<numSubTimers; ++i) {
   totMem  += timeData[i].memory;
   totTime += timeData[i].time;
 }

// get global total
#ifdef DISTRIBUTED
 double totMemory = (double) totMem;
 if(structCom) totMemory = structCom->globalSum(totMemory);
 totMem = (long) totMemory;
 if(structCom) totTime = structCom->globalSum(totTime);
 if(structCom) nSub = structCom->globalSum(nSub);
#endif 

 TimeData ret;
 ret.memory = totMem  / nSub;
 ret.time   = totTime / nSub;
 return ret;

}

TimeData
DistTimer::getTot()
{
 double totMemory = (double) overAll.memory;
 double totTime   = overAll.time;

#ifdef DISTRIBUTED
 if(structCom) totMemory = structCom->globalSum(totMemory);
 if(structCom) totTime   = structCom->globalMax(overAll.time);
#endif 

 TimeData ret;
 ret.memory = (long) totMemory; 
 ret.time   = totTime;
 return ret;

}

