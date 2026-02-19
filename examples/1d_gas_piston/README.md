# Full Genetic Optimization of 1-D piston subjected to repeated shock loading

## Piston

A simple spring-mass system solver is provided in this example repository to serve as the structural solver. To compile it, simply run

```sh
cd piston
make
```

The fluid solver used to capture the complex shock dynamics and fluid-structure interfactions for this problem is set to M2C.

## Dakota Setup

Asynchronous evaluations are spwaned using `Dakota`'s `fork` application interface, with the required setup specified in `dakota.in` input file. The `fork` interface requires an `analysis_driver` that reads the provided desgin parameters, performs the neccessary evaluations, and outputs the response functions. In `SOFICS`, the `driver.sh` bash script located in your `build` directory serves as the `analysis_dirver`. This script requires user-defined setup details, including input files for `Gmsh`, `M2C`, and `Aero-S`, as well as resource specifications for each evaluation. The setup details can be specified in a conguration file, where as the finite element mesh setup can be specified in a custom pre-processor script. Collectively, these scripts can be passed to the driver through command line arguments, like the ones shown in `dakota.in`. 

## Evaluation

To launch the simulation use `Dakota's` command line interface, i.e.,

```sh
dakota -i dakota.in -o dakota.log -w dakota.rst
```

Here, `dakota.in` is our input file, `dakota.log` is the log file to which Dakota's output will be written, and `dakota.rst` is a restart file.

## Cluster Evaluation

The `SLURM` scheduller is employed to launch the `Dakota` process on Virginia Tech's `Tinkercliffs` compute cluster. An example `SLURM` configuration can be found in the `run.sh` file. Update the following lines to match your preference and account details:

```sh
#SBATCH --job-name=dakota           # Job name
#SBATCH --partition=normal_q        # Partition or queue name
#SBATCH --account=m2clab            # Cluster account
```

The script uses the `dakota` command to call your `Dakota` installation, so ensure `Dakota` is properly installed before submitting a job. Follow the installation instructions available on the official [Dakota repository](https://github.com/snl-dakota/dakota?tab=coc-ov-file).

***Note:*** Ensure that sufficient compute nodes are allocated to the job. In this demonstration, each simulation requires 64 computational cores (CPUs). Therefore, the total number of cores needed will be `64 × evaluation concurrency`. Since each node on `Tinkercliffs` consists of 128 CPUs, you should adjust your resource allocation accordingly. Update the following line in `run.sh` to specify your resource requirements:

```sh
#SBATCH --nodes=4                   # Number of nodes
#SBATCH --ntasks-per-node=128       # Number of tasks per node
```

To submit a job on the compute cluster, use:

```sh
sbatch run.sh
```

To verify that the job was successfully submitted, run:

```sh
squeue | grep "your-user-id"
```

This command will display a list of jobs currently running under your user ID on the cluster.

## Results

