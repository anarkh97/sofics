ANDROID_NDK_ROOT := /home/avery/Android/android-ndk-r8d
NDK_TOOLCHAIN := arm-linux-androideabi-4.8.0
APP_MODULES = FemLib AcmeLib CommLib CorotationalLib DecLib DriverLib ElementLib FetiLib GNU-getoptLib HelmAxiLib HeteroLib LinpackLib MaterialLib MathLib MortarLib ParalLib ParserLib ProblemsLib RomLib SfemLib SolversLib ThreadsLib TimersLib UtilsLib
APP_STL := gnustl_static
APP_CFLAGS := -O3 -DNDEBUG
NDK_DEBUG := 0
