#ifndef DIST_TIMER_H_
#define DIST_TIMER_H_


struct TimeData {
  long memory;
  double time;
  TimeData() { memory = 0; time = 0.0; }
};

struct TimeStats {
  TimeData min;
  TimeData max;
  TimeData avg;
  TimeData tot;
  TimeStats(const TimeData &minT, const TimeData &maxT, 
            const TimeData &avgT, const TimeData &totT) : 
            min(minT), max(maxT), avg(avgT), tot(totT) { }
};

class DistTimer {
   TimeData overAll;
   TimeData *timeData;
   int numSubTimers;
  public:
    DistTimer(int _numSubTimers);
    ~DistTimer() { if(timeData) delete [] timeData; }
    void print();
    void printOverAll();

    const TimeData getMin();
    TimeData getMax();
    TimeData getAvg();
    TimeData getTot();
    TimeStats getStats();
    TimeData* getOverAll() { return &overAll; }

    int getNumSubTimers() { return numSubTimers; }

    void addTo(int iSub, long mem, double t) {  
       timeData[iSub].memory += mem; 
       timeData[iSub].time   += t; }
    void addOverAll(long mem, double t)
         { overAll.memory += mem; overAll.time += t; }
};

#endif
