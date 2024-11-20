#include <chrono>

#include <sys/time.h>
#include <Timers.d/GetTime.h>

using namespace std::chrono;
double getTime() {
	static auto t0 = high_resolution_clock::now();
	std::chrono::duration<double> elapsed_seconds = (high_resolution_clock::now()-t0);
	return elapsed_seconds.count()*1000;
}