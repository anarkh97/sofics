//
// Created by Michel Lesoinne on 2/22/18.
//

#ifndef FEM_AUTOTIMER_H
#define FEM_AUTOTIMER_H

#include <chrono>
#include <atomic>
#include <limits>

/** \brief Automatic timer system.
 * \details Use object of this type to easily time any desired operation.
 * By default, the timer start when it is constructed and stops at destruction. It keeps track of total time,
 * minimum time, maximum time and how many were destroyed.
 * Use a different id for each timer. To help with this point,  the _tid literal operator has been defined. Example:
 * AutoTimer<"computation"_tid> timer;
 * */
template <uint64_t id>
class AutoTimer {
public:
	struct Data {
		std::atomic<uint64_t> count{0}; //!< \brief Number of time a timer of this type was created.
		std::atomic<uint64_t> totalTime{0}; //!< \brief Cummulative time of all timers of this type.
		std::atomic<uint64_t> minTime{std::numeric_limits<uint64_t>::max()}; //!< \brief Minimum time for any timer of this type.
		std::atomic<uint64_t> maxTime{0}; //!< \brief Maximum time for any time of this type.
	};

	/** \brief Constructor for the timer.
	 *
	 * @param startRunning Whether the timer should start at construction.
	 */
	AutoTimer(bool startRunning = true) : isRunning(startRunning) {}

	~AutoTimer() {
		stop();
	}

	/** \brief Stop the timer.
	 *
	 * @param count Increment for number of call counter.
	 */
	void stop(int count = 1) {
		auto t1 = std::chrono::high_resolution_clock::now();
		if(!isRunning)
			return;
		isRunning = false;
		std::chrono::duration<uint64_t,std::nano> duration = t1-t0;
		Data &data = getModifiableData();
		data.count.fetch_add(1);
		data.totalTime.fetch_add(duration.count());
		uint64_t oldMin = data.minTime.load();
		while(oldMin > duration.count() && !data.minTime.compare_exchange_weak(oldMin, duration.count()))
			;
		uint64_t oldMax = data.maxTime.load();
		while(oldMax < duration.count() && !data.maxTime.compare_exchange_weak(oldMax, duration.count()))
			;
	}

	/** \brief Start the timer. */
	void start() {
		t0 = std::chrono::high_resolution_clock::now();
		isRunning = true;
	}

	static const Data &getData() {
		return getModifiableData();
	}
private:
	std::chrono::high_resolution_clock::time_point t0 = std::chrono::high_resolution_clock::now();
	bool isRunning = true;
	
	static Data &getModifiableData() {
		static Data data;
		return data;
	}
};

/** \brief Hash operator on string literals. Used for creating unique template types for easy to remember names. */
constexpr uint64_t operator "" _hash(const char *s, unsigned long length) {
	uint64_t result=0;
	for(size_t i = 0;i < length; ++i)
		result = result*((uint64_t )0x1134567) + ((uint64_t) s[i]);
	return result;
}


#endif //FEM_AUTOTIMER_H
