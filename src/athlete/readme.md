# Athlete project
- todo(lisca): update the description of the project

## 1. System dependencies

### 1.1 Ubuntu 22.04

todo(lisca): if you specify explicitly the compiler, can you use the same installation instructions for both Ubuntu 20.04 and Ubuntu 22.04?

### 1.2 Nvidia driver 535

Install the Nvidia driver:

```
sudo apt install nvidia-driver-535
```

<span style="color:red">
&rarr; Reboot your computer! Otherwise the tensorflow and pytorch will fail.
</span>

### 1.4 MuJoCo system dependencies.

Install MuJoCo's system dependencies.

```
sudo apt install \
    libx11-dev \
    libglew-dev \
    patchelf
```

### 1.4 OpenSim's system dependencies

OpenSim or BTK might not compile with the C++11 compiler. Therefore, we will ensure that for compilation, we will use an older compiler gcc-9. Install OpenSim's system dependencies:

#### 1.4.1 GCC 9.5.0 compiler 

```
sudo apt install g++-9 gcc-9
```

todo (lisca): Which C++ standard does GCC 9.5.0 compiler support?

todo (lisca): explain how to install multiple compilers.

https://www.linuxcapable.com/how-to-install-gcc-compiler-on-ubuntu-linux/


#### 1.4.2 OpenJDK 8

```
sudo apt install openjdk-8-jdk
```

#### 1.4.3 Matlab gtk-cambera

```
todo(lisca): update the installation instructions
```

## 2. Worspace layout

2.1 Export the `WS_BIOMEC` environment variable:

```
export WS_BIOMEC=$HOME/biomec
```

2.2 Create the following directory tree into your home directory:

```
cd $HOME && \
mkdir -p $WS_BIOMEC/build/btk/ && \
mkdir -p $WS_BIOMEC/build/opensim/core/ && \
mkdir -p $WS_BIOMEC/build/opensim/dependencies/ && \
mkdir -p $WS_BIOMEC/build/opensim/gui && \
mkdir -p $WS_BIOMEC/build/mujoco/ && \
mkdir -p $WS_BIOMEC/install/btk/ && \
mkdir -p $WS_BIOMEC/install/matlab/ && \
mkdir -p $WS_BIOMEC/install/mujoco/mujoco210/ && \
mkdir -p $WS_BIOMEC/install/opensim/core/ && \
mkdir -p $WS_BIOMEC/install/opensim/dependencies/ && \
mkdir -p $WS_BIOMEC/install/opensim/gui/ && \
mkdir -p $WS_BIOMEC/src/anaconda/ && \
mkdir -p $WS_BIOMEC/src/matlab/ && \
mkdir -p $WS_BIOMEC/src/netbeans/ && \
mkdir -p $WS_BIOMEC/data/
```

### 3 Matlab _R2022a_

&rarr; Note: _R2022_ because the package `xml_toolbox`, on which the BTK depends, has issues with Matlab versions newer than _R2022a_.

3.1 Download the Matlab instalation kit from [here](https://www.mathworks.com/downloads/) into the `WS_BIOMEC/src/matlab/` directory.

3.2 Unzip the downloaded kit and begin the installation:

```
cd $WS_BIOMEC/src/matlab/ && \
unzip matlab_R2022a_glnxa64.zip -d r2022a && \
cd r2022a && \
bash install
```

Enter your credentials (email and password) for the Matworks website. Click on the `Sign in` button, and accept the license agreement, hit the `Next` button twice.

In the step `DESTINATION`, for the `Select destination folder`, enter the following string:

```
/home/<use name>/biomec/install/matlab/r2022a/
```

and replace the part <user name> with your username. Click the `Next` button.

Select the additional packages which we need to install:

&rarr; [DSP System Toolbox](https://de.mathworks.com/products/dsp-system.html)

&rarr; [Signal Processing Toolbox](https://de.mathworks.com/products/signal.html)

&rarr; [Statistics and Machine Learning Toolbox](https://www.mathworks.com/products/statistics.html)

Click the `Next` button.

Ensure that the check box for `Create symbolic links to MATLAB scripts in:` is unchecked!

Click `Next`, and afterward `Begin Install`.

Your installation will proceed.

Feel free to grab yourself a cup of coffee, because depending on your internet speed, it might take a while until Matlab downloads and installs all packages which it needs.

3.3 Test your installation:

Open a new shell (terminal) and run the follwing commands.

Run the step 2.1, which exports the `WS_BIOMEC` system variable.

```
$WS_BIOMEC/install/matlab/r2022a/bin/matlab
```

If  your Matlab r2022a starts, then your installation was successful. If Matlab does not start, then please report this to the developers of this repository.

3.4 Remove the `libstdc++` which came with Matlab, because it crashes Matlab, because it its different than the `libstdc++` which came with Ubuntu 22.04.

```
rm $WS_BIOMEC/install/matlab/r2022a/sys/os/glnxa64/libstdc++.so.6
```


## 4. Anaconda

The entire software stack will be contained into conda environment, isolated from your system packages. Only a few debian packages will be installed into your system.

4.1 Download Anaconda.

```
cd $WS_BIOMEC/src/anaconda/ && \
wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh
```

4.2 Install Anaconda.

```
bash $WS_BIOMEC/src/anaconda/Anaconda3-2024.06-1-Linux-x86_64.sh -p $WS_BIOMEC/install/a202406/
```

Confirm the installation process:

```
Welcome to Anaconda3 2024.06-1

In order to continue the installation process, please review the license
agreement.
Please, press ENTER to continue
>>> < just press the ENTER key >
```

Accept the license agreement:

```
...
...
Do you accept the license terms? [yes|no]
```

Type in `yes` and press enter.

```
[no] >>> yes
```

Confirm the location:
```
...

...

Anaconda3 will now be installed into this location:
/home/<your use name>/biomec/install/a202406/

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/<your user name>/biomec/install/a202406/] 
```

Press enter.

```
>>> < enter key >
```

&rarr; Feel free to grab yourself a coffee. Probably it will take Anaconda more tha 2 minutes (depending on your internet speed) to to download and install all default packages which it needs.

When the installer asks:

```
...

...

Do you wish the installer to initialize Anaconda3
by running conda init? [yes|no]
```

Type in `no` and press enter, because in the future you will initialize anaconda manually, otherwise anaconda will be activated automatically each time you start a new temrinal. We definitely do not want that!

You are done with installing Anaconda3.

## 5. `athlete_py38` conda environment

5.1 Clone the `athlete` (this repository) repository inside the `src` directory.

```
cd $WS_BIOMEC/src/ && \
git clone git@gitlab.lrz.de:glisca/athlete.git
```

5.2 Symbolicaly link `athlete`'s local directory `data/local/data/` to the `$HOME/biomec/data/` directory. In this way the logs of the RL algorithms will be saved into the `$HOME/biomec/data/` directory, and will not polute the `athlete` repository.

```
cd $WS_BIOMEC/src/athlete/ && \
mkdir -p data/local/ && \
ln -s $WS_BIOMEC/data data/local/data
```

5.3 Create the `athlete_py38` conda environment (and automatically install the required dependencies).

&rarr; Run the following command and feel free to grab yourself a coffee. Probably it will take anaconda more tha 5 minutes (depending on your internet speed) to to download and install all packages which the about to be installed `athlete_py38` conda evironment needs.

```
$WS_BIOMEC/install/a202406/bin/conda env create -f $WS_BIOMEC/src/athlete/athlete/config/conda/athlete_py38.yml
```

A message similar to the following:
```
...
done
#
# To activate this environment, use
#
#     $ conda activate athlete_py38
#
#     $ conda deactivate

Retrieving notices: ...working... done
```
confirms that your new `athlete_py38` conda environment was successfully installed.

5.6 Activate manually the `athlete_py38` conda environment.

```
source ~/biomec/src/athlete/athlete/bin/activate_athlete_py38.sh
```

## 6. [MuJoCo](https://mujoco.org/)

We will use an older version of [MuJoCo](https://mujoco.org/), namely the one on which `garage` depends.

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 6. assume that you have the `athlete_py38` conda environment activated.
</span>

6.1 Install [MuJoCo](https://mujoco.org/).

```
cd $WS_BIOMEC/install/mujoco/ && \
wget https://mujoco.org/download/mujoco210-linux-x86_64.tar.gz && \
tar -xvf mujoco210-linux-x86_64.tar.gz && \
cd mujoco210/bin && \
wget https://www.roboti.us/file/mjkey.txt
```

6.2 Test [MuJoCo](https://mujoco.org/)'s installation.

```
cd $WS_BIOMEC/install/mujoco/mujoco210/bin && \
./simulate ../model/humanoid.xml
```

6.3 Finalize [mujoco-py](https://github.com/openai/mujoco-py)'s installation by cytonizing MuJoCo's Python bindings:

[mujoco-py](https://github.com/openai/mujoco-py) has been installed during the creation of the `athlete_py38` conda environment. By typing the following commands, it will be compiled and configured, such that it can be used together by `garage`.

```
cd $WS_BIOMEC/src/athlete && \
python athlete/run/test_mujoco_installation.py
```

No crash during cythonization and an output similar to:

```
[0.  0.  1.4 1.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
 0.  0.  0.  0.  0.  0.  0.  0.  0.  0. ]
[-1.12164337e-05  7.29847036e-22  1.39975300e+00  9.99999999e-01
  1.80085466e-21  4.45933954e-05 -2.70143345e-20  1.30126513e-19
 -4.63561234e-05 -1.88020744e-20 -2.24492958e-06  4.79357124e-05
 -6.38208396e-04 -1.61130312e-03 -1.37554006e-03  5.54173825e-05
 -2.24492958e-06  4.79357124e-05 -6.38208396e-04 -1.61130312e-03
 -1.37554006e-03 -5.54173825e-05 -5.73572648e-05  7.63833991e-05
 -2.12765194e-05  5.73572648e-05 -7.63833991e-05 -2.12765194e-05]
```

are good signs! :). Your 1st forward simulation in MuJoCo is successful.

```
quit()
```

<span style="color:red">
&rarr; 6.4 (! ! ! skip this step) Building from source
</span>

```
sudo apt install libxrandr-dev libxinerama-dev libxcursor-dev 
```

```
cmake --build . -j$(nproc)
```

<span style="color:red">
&rarr;
6.5 (! ! ! skip this step) Install from binaries and build the Python bidndings from source
</span>

Just follow MuJoCo's [Python / Building from source](https://mujoco.readthedocs.io/en/latest/python.html#building-from-source) instructions.

! MuJoCo 3.2.4 dropped the support for Python 3.8.

## 7. [garage](https://garage.readthedocs.io/)

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 7. assume that you have the `athlete_py38` conda environment activated.
</span>

todo(lisca): insert the commands to install manually the `tensorflow` and `pytorch` libraries!

7.0 Test the tensorflow and pytorch of the `athlete_py38` conda environment:

Start python.

```
cd $WS_BIOMEC/src/athlete && \
python athlete/run/test_tensorflow_and_pytorch.py
```

If the creation of the `athlete_py38` conda environment was successfull then the last two lines of the output should look like:

```
tensorflow report: True
pytorch    report: True
```

If you got here successfully then your Machine Learning software stack is ready.

7.0 Install the gym==0.19.0 (the version garage requires) according to these [instructions](https://stackoverflow.com/questions/77124879/pip-extras-require-must-be-a-dictionary-whose-values-are-strings-or-lists-of/77205046). The instructions explain the need to downgrade `setuptools`, `pip` and `wheel`.

```
pip install setuptools==65.5.0 pip==21 wheel==0.38.0 && \
pip install gym==0.19.0
```

7.1 Clone the garage's repository and install it:

```
cd $WS_BIOMEC/src/ && \
git clone https://github.com/rlworkgroup/garage.git && \
cd garage/ && \
pip install --upgrade --upgrade-strategy only-if-needed -e .
```

7.2 Test [garage](https://garage.readthedocs.io/)'s installation by running a few RL algorithms:

BC (Behavioral Cloning)

```
cd $WS_BIOMEC/src/garage && \
python src/garage/examples/torch/bc_point_deterministic_policy.py
```

In the console you should see printed statistics about the progress of the training.

To stop the training, press `Ctrl+c`. 

DDPG

```
cd $WS_BIOMEC/src/garage && \
python src/garage/examples/tf/ddpg_pendulum.py
```

TD3

```
cd $WS_BIOMEC/src/garage && \
python src/garage/examples/tf/td3_pendulum.py
```

PPO

```
cd $WS_BIOMEC/src/garage && \
python src/garage/examples/tf/ppo_pendulum.py
```

SAC

```
cd $WS_BIOMEC/src/garage && \
python src/garage/examples/torch/sac_half_cheetah_batch.py
```

7.3 Visualize the progress of the trainings in `tensorboard`

Open a new terminal and activate the `athlete_py38` conda environment.

```
export WS_BIOMEC=$HOME/biomec && \
cd $WS_BIOMEC/src/garage && \
source $WS_BIOMEC/src/athlete/athlete/bin/activate_athlete_py38.sh
```

Start the tensorboard:

```
tensorboard --logdir data/local/experiment/
```

In your browser, open the following address:

```
http://localhost:6006/
```

7.4 [OpenAI's Spinning Up in Deep RL!](https://spinningup.openai.com/en/latest/#welcome-to-spinning-up-in-deep-rl)

Read this introduction in the field of Deep Reinforcement Learning which summarizes the most important deep reinforcement learning algorithms.


## 8. [OpenSim core](https://github.com/opensim-org/opensim-core)

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 8. assume that you have the `athlete_py38` conda environment activated.
</span>

8.1 Install OpenSim 4.5.1 system dependencies:

```
sudo apt install \
    build-essential \
    gcc-9 \
    g++-9 \
    libtool \
    autoconf \
    pkg-config \
    gfortran \
    ninja-build \
    liblapack-dev \
    freeglut3-dev \
    libxi-dev \
    libxmu-dev
```

8.2 Clone the OpenSim repository into the 'src' directory and checkout the specific tag:

```
cd $WS_BIOMEC/src/ && \
git clone https://github.com/opensim-org/opensim-core.git && \
cd opensim-core && \
git checkout 4.5.1
```

8.3 Build and (automatically) install OpenSim's dependencies.

Configure the build:

```
cd $WS_BIOMEC/build/opensim/dependencies/ && \
cmake $WS_BIOMEC/src/opensim-core/dependencies/ \
    -DCMAKE_INSTALL_PREFIX=$WS_BIOMEC/install/opensim/dependencies/ \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++-9 \
    -DCMAKE_CXX_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_CXX_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DCMAKE_C_COMPILER=/usr/bin/gcc-9 \
    -DCMAKE_C_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_C_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DOPENSIM_WITH_CASADI=ON \
    -DOPENSIM_WITH_TROPTER=ON \
    -DSUPERBUILD_ezc3d=ON
```

Build and install the binaries:

```
cmake --build . -j$(nproc)
```

A message like the following:

```
...

[ 97%] Completed 'adolc'
[100%] Built target adolc
```

It's a good indicator that the compilation OpenSim's dependencies was successful.

8.4 Build OpenSim's binaries

We will generate the Matlab bindings for OpenSim too. Therefore, ensure that the `$BIOMEC_MATLAB_ROOT` bash variable points to the root directory where your Matlab is installed:

```
echo $BIOMEC_MATLAB_ROOT
```

Configure the build:

```
cd $WS_BIOMEC/build/opensim/core/ && \
cmake $WS_BIOMEC/src/opensim-core/ \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX=$WS_BIOMEC/install/opensim/core/ \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++-9 \
    -DCMAKE_CXX_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_CXX_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DCMAKE_C_COMPILER=/usr/bin/gcc-9 \
    -DCMAKE_C_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_C_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DBUILD_TESTING=OFF \
    -DBUILD_JAVA_WRAPPING=ON \
    -DBUILD_PYTHON_WRAPPING=ON \
    -DOPENSIM_C3D_PARSER=ezc3d \
    -DOPENSIM_DEPENDENCIES_DIR=$WS_BIOMEC/install/opensim/dependencies/ \
    -DOPENSIM_INSTALL_UNIX_FHS=OFF \
    -DOPENSIM_PYTHON_STANDALONE=ON \
    -DMatlab_ROOT_DIR=$BIOMEC_MATLAB_ROOT
```

Compile the configured build:

```
cmake --build . -j$(nproc)
```

<span style="color:red">
&rarr; repeat the previous command until no "collored" message is printed in the console.
</span>

Without a complete comilation of all its binaries OpenSim will probably malfunction.

A message similar to the following:

```
...
[ 96%] Built target PythonBindings
[100%] Built target osimJavaJNI
[100%] Built target JavaBindings
```
is a good indication that the compilation of the OpenSim binaries (executables and libraries) was successful.

8.5 Install OpenSim's binaries:

```
cmake --install .
```

8.6 Install the Python bindings of OpenSim into the `athlete_py38` conda environment.

```
cd $WS_BIOMEC/install/opensim/core/sdk/Python/ && \
pip install .
```

8.7 Test OpenSim's installation into the `athlete_py38` conda environment.

```
cd $WS_BIOMEC/src/athlete && \
python athlete/run/test_opensim_installation.py
```

## 9. OpenSim [GUI](https://github.com/opensim-org/opensim-gui) following these [build instructions](https://github.com/opensim-org/opensim-gui/blob/main/scripts/build/opensim-gui-linux-build-script.sh):

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 9. assume that you have the `athlete_py38` conda environment activated.
</span>

9.1 Ensure that the following Debian packages are installed:

```
sudo apt install \
    build-essential \
    autotools-dev \
    autoconf \
    pkg-config \
    automake \
    liblapack-dev \
    freeglut3-dev \
    libxi-dev \
    libxmu-dev \
    doxygen \
    python3 \
    libpcre3 \
    libpcre3-dev \
    libpcre2-dev \
    byacc \
    gfortran \
    libtool \
    libssl-dev \
    libffi-dev \
    ninja-build \
    patchelf
```

9.2 Download and install NetBeans 12.3.

```
cd $WS_BIOMEC/src/netbeans/ && \
wget -nc -q --show-progress https://archive.apache.org/dist/netbeans/netbeans/12.3/Apache-NetBeans-12.3-bin-linux-x64.sh && \
bash Apache-NetBeans-12.3-bin-linux-x64.sh
```

When the GUI of Netbeans, click the `Next` button. Enter the following path:

```
/home/<use name>/biomec/install/netbeans-12.3
```

into the `Install the Apache NetBeans IDE to:` variable and click the `Next`, `Install` and `Finish` buttons, as they appear.

9.3 Clone the `opensim-gui` repository.

&rarr; Run the following command and feel free to grab yourself a coffee. Probably it will take git more than 5 minutes (depending on your internet speed) to to download the OpenSim GUI repository from GitHub.

```
cd $WS_BIOMEC/src/ && \
git clone https://github.com/opensim-org/opensim-gui.git && \
cd opensim-gui && \
git checkout 4.5 && \
git submodule update --init --recursive -- opensim-models opensim-visualizer Gui/opensim/threejs 
```

9.4 Build OpenSim's GUI

Configure the build:

```
cd $WS_BIOMEC/build/opensim/gui && \
cmake $WS_BIOMEC/src/opensim-gui/ \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_INSTALL_PREFIX=$WS_BIOMEC/install/opensim/gui/ \
    -DCMAKE_PREFIX_PATH=$WS_BIOMEC/src/opensim-core/ \
    -DAnt_EXECUTABLE=$WS_BIOMEC/install/netbeans-12.3/netbeans/extide/ant/bin/ant \
    -DANT_ARGS="-Dnbplatform.default.netbeans.dest.dir=$WS_BIOMEC/install/netbeans-12.3/netbeans;-Dnbplatform.default.harness.dir=$WS_BIOMEC/install/netbeans-12.3/netbeans/harness"
```

Compile the configured build:

```
make CopyOpenSimCore -j$(npro) && \
make PrepareInstaller -j$(nproc)
```

9.5 Install the fresh build

```
cd $WS_BIOMEC/src/opensim-gui/Gui/opensim/dist/installer/opensim
```

Hack the INSTALL script that your OS is Ubuntu 18.04. Open the directory containing the `INSTALL` script.

```
nautilus .
```

Open the `INSTALL` script with a text editor, and replace the lines 87 - 88 with the following ones:

```
    OS="Ubuntu" # $NAME
    VER="18.04" # $VERSION_ID
```

Now, the script will thik that your Ubuntu is Ubuntu 18.04.

In the shell, run the installer and let the script use sudo, even the install prefix should not need sudo. Aparently the script installs some additional packages essential for OpenSim's gui.

```
bash INSTALL --prefix=$WS_BIOMEC/install/opensim/gui/
```

Thest your installation:

```
opensim_gui
```

This is an alias to the binary of the OpenSim GUI. The `athlete/bin/activate_athlete_py38.sh` script is responsible for exporting it correctly.

## 10. [BTK](http://biomechanical-toolkit.github.io/)

We use this library in order to process the Vicon data.

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 11. assume that you have the `athlete_py38` conda environment activated.
</span>

10.1 Clone the repository

```
cd $WS_BIOMEC/src/ && \
git clone https://github.com/Biomechanical-ToolKit/BTKCore.git btkcore
```

Configure the build:

```
cd $WS_BIOMEC/build/btk/ && \
ccmake $WS_BIOMEC/src/btkcore/ \
    -DCMAKE_INSTALL_PREFIX=$WS_BIOMEC/install/btk/ \
    -DCMAKE_CXX_COMPILER=/usr/bin/g++-9 \
    -DCMAKE_CXX_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_CXX_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DCMAKE_C_COMPILER=/usr/bin/gcc-9 \
    -DCMAKE_C_COMPILER_AR=/usr/bin/gcc-ar-9 \
    -DCMAKE_C_COMPILER_RANLIB=/usr/bin/gcc-ranlib-9 \
    -DBTK_WRAP_MATLAB=ON \
    -DMATLAB_ROOT=$BIOMEC_MATLAB_ROOT
```

Press the key `c` to configure. Check if there are error messages. Warning messages are OK.

Press the key `e` to return to the list of the variables. The key `g` is not yet available.

Press (again) the key `c` to configure.

Press the key `e` to return to the list of the variables.

Now the key `g` became available.

Press the key `g` to generate the config files for the build process.

10.2 Build

```
cmake --build . -j$(nproc)
```

10.3 Install

```
cmake --install .
```

&rarr; If you want to build the Python bindings then you must manually set the following variables too:

```
-DNUMPY_INCLUDE_DIR=<path to the include directory of the numpy>
-DNUMPY_VERSION=<version of the numpy>
-DBTK_WRAP_PYTHON=ON

```

You can find them by runnin the following Python code snipped into your conda envirohment:

```
import numpy as np

print(np.get_include())
print(np.version.version)
```

## 11. Make OpenSim's installation visible to Matlab

<span style="color:red">
&rarr; Make sure that your `athlete_py38` conda environment is activated with the command from the section 5.6.
</span>

<span style="color:red">
&rarr; All substeps of step 11. assume that you have the `athlete_py38` conda environment activated.
</span>

<span style="color:red">
&rarr; If the `athlete_py38` conda environment is not activated, then (OpenSim Matlab Tools)[] will not find the `opensim-cmd run-tool` binary!
</span>

11.1 Start Matlab (using the alias) without the GUI:

```
cd $WS_BIOMEC/install/opensim/core/Resources/Code/Matlab/ &&\
matlabr2022a -nodisplay
```

11.2 Make OpenSim visible from Matlab.

From within Matlab, run the configureOpenSim script, and exit:
```
configureOpenSim();
quit();
```

You should recognize an output similar to:

```
>> configureOpenSim();
-- Added /home/<your user name>/biomec/install/opensim/core/sdk/Java/org-opensim-modeling.jar to /home/lisca/.matlab/R2022a/javaclasspath.txt.

-- Added /home/<your user name>/biomec/install/opensim/core/sdk/lib to /home/lisca/.matlab/R2022a/javalibrarypath.txt.

-- Added /home/<your user name>/biomec/install/opensim/core/Resources/Code/Matlab/Utilities to the MATLAB path.

>>quit();
```

Hit enter, and Matlab will exit.

11.3 Test the installation, by starting Matlab (using the alias) without the GUI, and running the dummy code snipped:

```
import org.opensim.modeling.*;
model = Model();
```

is a good sign that the configuration of Matlab to use OpenSim was successful.

Stop Matlab and start it from from the same terminal again:

```
matlabr2022a -nodisplay
```

Test the installation with the lines 5 - 6 of the `$WS_BIOMEC/install/opensim/core/Resources/Code/Matlab/configureOpenSim.m` which you previously run to configure / make visible OpenSim into Matlab:

```
import org.opensim.modeling.*;
model = Model();
```

"No news are good news." If no error is returned then the configuration of Matlab to use OpenSim was successful.

At this point the installation is done. The rest of the sections are still present only for the purpose to be cleaned up and refactored.

## 13 Issues

### 13.1 `libstdc++`

Problem: Matlab crashes:

```
/home/lisca/biomec/install/matlab/r2022a/bin/glnxa64/cefhelper: /home/lisca/biomec/install/matlab/r2022a/bin/glnxa64/../../sys/os/glnxa64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found (required by /home/lisca/biomec/install/a202406/envs/athlete_py38/lib/././libicuuc.so.73)
```

Solution ideas:

1. simply delete the libraries that come with Matlab, as recommended [here](https://de.mathworks.com/matlabcentral/answers/1907290-how-to-manually-select-the-libstdc-library-to-use-to-resolve-a-version-glibcxx_-not-found)

2. Install the version of the libstdcxx-ng which is the same as the version of the `libstdc++6` / `libstd++-9-dev` / `libstdc++-11-dev` Ubuntu22.04.

### 13.2 `ncurses`

Problem:

```
tig: /home/lisca/biomec/install/a202406/envs/athlete_py38/lib/libtinfo.so.6: no version information available (required by tig)
tig: /home/lisca/biomec/install/a202406/envs/athlete_py38/lib/libtinfo.so.6: no version information available (required by tig)
tig: /home/lisca/biomec/install/a202406/envs/athlete_py38/lib/libncursesw.so.6: no version information available (required by tig)
```

Solution ideas:

1. install `ncurses==6.3` from conda, because Ubuntu22.04 has `libncurses6:amd64 6.3-2ubuntu0.1` installed.

## 14 MyoSuite

```
conda install -c conda-forge libstdcxx-ng
```

## 15 Graveyard

## 9. Athlete: [OpenSim pipeline tools](https://simtk.org/projects/matlab_tools)

Make sure that your `athlete_py38` conda environment is activated, with the command from the section 4.4.

Start Matlab:

```
matlab -softwareopengl
```

Add the directory `walking_pipeline` and its subdirectories to Matlab's `Path`.

Add the directory `$HOME//biomec/sim/install/btk/` and its subdirectories to Matlab's `Path`.

Open the `matlab/buildLegModel.m` file and run it.

Open the `matlab/walking_pipeline.m` file and run it.

Choose the file `squad_side_0_static.c3d`.

Run the script again, and this time choose the file `gait_1x_0.c3d`.

## 10. Athlete: Moco

Make sure that your `athlete_py38` conda environment is activated with the command from the section 4.4.

Install the athlete package into the `athlete_py38` conda environment.

```
cd $HOME/biomec/sim/src/athlete/ && \
pip install .
```

Run OpenSim Moco on the default model.

```
python athlete/run/optimizer.py --optimize
```

Visualize the results.

```
python athlete/run/optimizer.py --simulate
```

## 11. Athlete: Reinforcement Learning

Start the training.

```
python athlete/run/manager.py
```

The manager will reinstall the `garage` and `athlete` packages and start the training. In the console you will see details about progress of the installation and training.

For visualizing the musculoskeletal model durig training, ask the owner of the repository for indications.

# Troubleshooting


## ?. 

<span style="color: orange">
&rarr; Some binaries (libraries or executable) might not be build during the 1st run of the previous command. </span> Therefore, run the last command a few times more until the building of all binnaries is successful.


<b><span style="color: red">
&rarr; Only if the building still fails! </span></b> then find out the number of processing units which your CPU has:


```
nproc
```

Delete the old build:

```
cd $HOME/biomec/sim/build/core/ && \
rm -rf *
```

Configure a new build as in the 1s step of this section.

Build again the new configuration by using only __half__ of the processing units which CPU has (half of the output of the `nproc` command).

```
cmake --build . <replace the brakets with half of the number returned by the `nproc` command>
```


## 1. The training crashes because of the `ray` library.

If the training crashes due to the default version of ray, which garage installs, then install the version [ray-2.0.0.dev](https://docs.ray.io/) compiled with Python 3.8. For more than a year, this version run stable.

Make sure that your `athlete_py38` conda environment is activated with the command from the section 4.4.

```
cd /tmp/ && \
wget https://s3-us-west-2.amazonaws.com/ray-wheels/latest/ray-2.0.0.dev0-cp38-cp38-manylinux2014_x86_64.whl && \
pip install ray-2.0.0.dev0-cp38-cp38-manylinux2014_x86_64.whl
```

## 2. Moco on models with locked joints

Moco can not handle locked joints. If you want to lock joints in your model, then you must convert them into WeldJoint.


# Legacy [mujoco-py](https://github.com/openai/mujoco-py) installation instructions

The following instructions are snippents of installation commands. They are not expected to run sequentially, rather remember of key steps and aspects of the installation process.

1. Install anaconda

```
wget https://repo.anaconda.com/archive/Anaconda3-2021.11-Linux-x86_64.sh
chmod +x Anaconda3-2021.11-Linux-x86_64.sh
./Anaconda3-2021.11-Linux-x86_64.sh
```

2. Download the Mujoco library from 

```
wget https://mujoco.org/download/mujoco210-linux-x86_64.tar.gz
mkdir /home/username/.mujoco
cd /home/username/.mujoco/
tar -xvf mujoco210-linux-x86_64.tar.gz
```

3. include these lines in  .bashrc file:

```
export LD_LIBRARY_PATH=/home/user_name/.mujoco/mujoco210/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/lib/nvidia
export PATH="$LD_LIBRARY_PATH:$PATH"
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libGLEW.so
```

4. Test that the library is installed by going into:

```
cd ~/.mujoco/mujoco210/bin
./simulate ../model/humanoid.xml
```

5. Install mujoco-py:
```
conda create --name mujoco_py python=3.8
conda activate mujoco_py
sudo apt update
sudo apt-get install patchelf
sudo apt-get install python3-dev build-essential libssl-dev libffi-dev libxml2-dev
sudo apt-get install libxslt1-dev zlib1g-dev libglew1.5 libglew-dev python3-pip

git clone https://github.com/openai/mujoco-py
cd mujoco-py
pip install -r requirements.txt
pip install -r requirements.dev.txt
pip install -e . --no-cache

```

6. reboot your machine

7. run these commands

```
conda activate mujoco_py
sudo apt install libosmesa6-dev libgl1-mesa-glx libglfw3
sudo ln -s /usr/lib/x86_64-linux-gnu/libGL.so.1 /usr/lib/x86_64-linux-gnu/libGL.so
python3 examples/setting_state.py
```
#
conda install -c conda-forge cudatoolkit=11.3.1 cudatoolkit-dev=11.3.1 cudnn=8.2.1

pip install cython==0.29.32 mujoco-py==2.1.2.14
pip install protobuf==3.20.*
pip install tensorboard==2.6.0 tensorboardx==2.4 tensorflow==2.6.0
pip install torch==1.7.1 torchvision==0.8.2
