from conans import ConanFile, tools, AutoToolsBuildEnvironment
import os


class OpenMPIConan(ConanFile):
    name = "OpenMPI"
    lower_name = "openmpi"
    version = "3.0.0"
    url = "https://github.com/bincrafters/conan-openmpi"
    description = "A High Performance Message Passing Library"
    license = "https://www.open-mpi.org/community/license.php"
    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False],
               "fortran": ['yes', 'mpifh', 'usempi', 'usempi80', 'no']}
    default_options = "shared=False", "fortran=no"

    #def requirements(self):
    #    self.requires.add("zlib/1.2.11@conan/stable")

    def system_requirements(self):
        if self.settings.os == "Linux":
            if tools.os_info.linux_distro == "ubuntu" or tools.os_info.linux_distro == "debian":
                installer = tools.SystemPackageTool()
                installer.install('openssh-client')

    def source(self):
        version_tokens = self.version.split('.')
        version_short = 'v%s.%s' % (version_tokens[0], version_tokens[1])
        source_url = "https://www.open-mpi.org/software/ompi"
        # https://www.open-mpi.org/software/ompi/v3.0/downloads/openmpi-3.0.0.tar.bz2
        download_url = "{0}/{1}/downloads/{2}-{3}.tar.bz2".format(source_url, version_short, self.lower_name, self.version)
        zip_name = os.path.basename(download_url)
        if(not os.path.isfile(zip_name)):
            tools.get(download_url)
        else:
            print("Found the source file in house")
        extracted_dir = self.lower_name + "-" + self.version
        os.rename(extracted_dir, "sources")

    def build(self):
        with tools.chdir("sources"):
            env_build = AutoToolsBuildEnvironment(self)
            args = ['--disable-wrapper-rpath',
                    'prefix=%s' % self.package_folder]
            if self.settings.build_type == 'Debug':
                args.append('--enable-debug')
            if self.options.shared:
                args.extend(['--enable-shared', '--disable-static'])
            else:
                args.extend(['--enable-static', '--disable-shared'])
            # args.append('--enable-mpi-fortran=%s' % str(self.options.fortran))
            # args.append('--with-zlib=%s' % self.deps_cpp_info['zlib'].rootpath)
            # args.append('--with-zlib-libdir=%s' % self.deps_cpp_info['zlib'].lib_paths[0])
            env_build.configure(args=args)
            env_build.make(args=['-j8'])
            env_build.make(args=['install', '-j8'])

    def package(self):
        with tools.chdir("sources"):
            self.copy(pattern="LICENSE")

    def package_info(self):
        self.cpp_info.libs = ['mpi', 'open-rte', 'open-pal']
        if self.settings.os == "Linux":
            self.cpp_info.libs.extend(['dl', 'pthread', 'rt', 'util'])
        self.env_info.MPI_HOME = self.package_folder
        self.env_info.MPI_BIN = os.path.join(self.package_folder, 'bin')
