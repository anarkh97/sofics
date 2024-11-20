from conans import ConanFile, CMake

class FetiLibConan(ConanFile) :
    settings = "os", "compiler", "build_type", "arch"
    # requires = "MPICH/3.2@cmsoft/stable"
    generators = "cmake"
    options = {"mpi": ["mpich", "openmpi"]}
    default_options = "mpi=openmpi"
    #requires = "OpenMPI/3.0.0@cmsoft/stable"

    def requirements(self):
        print ("Go with ", self.options.mpi)
        if self.options.mpi == "mpich":
            self.requires("MPICH/3.2@cmsoft/stable")
        else:
            self.requires("OpenMPI/3.0.0@cmsoft/stable")

    def imports(self):
        if self.options.mpi == "mpich":
            self.copy("hydra_pmi_proxy*", dst="bin", src="bin")  #for MPICH
            self.copy("mpiexec.hydra", dst="bin", src="bin") #for MPICH
        else:
            self.copy("opal_wrapper", dst="bin", src="bin")
        self.copy("mpic++", dst="bin", src="bin")
        self.copy("mpicc", dst="bin", src="bin")
        self.copy("mpicxx", dst="bin", src="bin")
        self.copy("orterun", dst="bin", src="bin")
        self.copy("mpirun", dst="bin", src="bin")
        self.copy("mpiexec", dst="bin", src="bin")
