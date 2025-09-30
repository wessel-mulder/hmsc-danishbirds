# finding python
system_python = "python3"
# system_python = "/Users/username/opt/anaconda3/envs/tf/bin/python3"
system2(system_python, "--version")

# set up hmsc
system2(system_python, "-m venv hmsc-venv")
python = file.path(getwd(), "hmsc-venv", "bin", "python")  # for Linux and macOS
# python = file.path(getwd(), "hmsc-venv", "Scripts", "python")  # for Windows
system2(python, "--version")

# install hmsc-hpc
package_path = "/projects/cmec/people/bhr597/projects/hmsc-danishbirds/hmsc-hpc-main"
system2(python, "-m pip install --upgrade pip")
system2(python, paste("-m pip install", shQuote(package_path)))

# Choose correct python by uncommenting correct line:
# python = "python3"  # default
python = file.path(getwd(), "hmsc-venv", "bin", "python")  # hmsc-venv for Linux and macOS
# python = file.path(getwd(), "hmsc-venv", "Scripts", "python")  # hmsc-venv for Windows

Sys.setenv(TF_CPP_MIN_LOG_LEVEL=3)  # reduce debug output from tensorflow
system2(python, "-c \"import tensorflow as tf; print(tf.constant(1))\"")
system2(python, "-c \"import hmsc\"")

