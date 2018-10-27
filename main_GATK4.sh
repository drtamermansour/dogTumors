# https://software.broadinstitute.org/gatk/documentation/quickstart

# 1. Requirements
module load Java/jdk1.8.0 ## java -version ## java version "1.8.0_102"
# 2. Get GATK
# 2A. download: https://software.broadinstitute.org/gatk/download/
wget https://github.com/broadinstitute/gatk/releases/download/4.0.9.0/gatk-4.0.9.0.zip
unzip gatk-4.0.9.0.zip
cd gatk-4.0.9.0
# 2B. Add the path to the .bashrc file "export PATH=$PATH:/mnt/home/mansourt/gatk-4.0.9.0/gatk"
source $HOME/.bashrc
# 2B. run the conda env (I have installed miniconda already & I have miniconda3/bin in the PATH)
# https://software.broadinstitute.org/gatk/documentation/article?id=12836
conda update conda  ## conda --version ## conda 4.5.11
conda env create -n gatk -f gatkcondaenv.yml
source activate gatk ## the test command in this documentation pagae failed to run! ## A USER ERROR has occurred: 'NeuralNetInference' is not a valid command.




