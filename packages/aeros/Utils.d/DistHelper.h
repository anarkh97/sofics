#ifndef DIST_HELPER_H_
#define DIST_HELPER_H_

#include <cstdio>
#include <cstdarg>
#include <set>

extern int salinasFlag;
extern int verboseFlag;

void filePrint(FILE *file, const char *format, ...);

bool imPrinting();

template <typename ...T>
void
filePrintOnce(FILE *file, const char *format, T && ...args)
{
	static std::set<const char *> seen;
	if(!salinasFlag || (salinasFlag && verboseFlag)) {
		if (seen.insert(format).second) {
#ifdef DISTRIBUTED
			if (imPrinting())
#endif
//#ifdef ANDROID
//  __android_log_print(ANDROID_LOG_ERROR, "FemLib", format, args);
//#else
				fprintf(file, format, args...);
//#endif
		}
	}
}
#endif
