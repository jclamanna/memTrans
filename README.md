# Membrane Transfer Classifier

## Description
Membrane Transfer Classifier is a program designed to process membrane and trace identification data from membrane transfer images that have been processed through ImageJ. The program classifies traces into the following categories:
- Docking Event
- Central Channel Abortive
- Cytoplasmic Fibril Abortive
- Successful
- Crossed Midline

Additionally, there is an option to plot all traces to view their trajectories through the membrane.

## Installation

To install the program:

1. Clone the repository.
2. Ensure you have Anaconda installed on your system.
3. Create a new conda environment using the provided `mem.yml` file with the command: 

`conda env create -f mem.yml`

4. Once the environment is set up, you are ready to use the program.

## Usage

To use the program:

1. Activate the conda environment with the command: `conda activate mem`
2. Run the program using the following command:

`python controller.py <Membrane File Path> <Trace File Path> <Pixel Size>`

For example:

`python controller.py data/membraneTestBed.csv data/tb_fit_results_2.csv 240`

This command will perform a basic run and generate an `out.csv` file containing all classified traces and some basic statistics about them.

### Optional Arguments

- `-p`: This argument will plot all traces.
- `-h`, `--help`: Use this argument to display help information.
