# Goodman High Throughput Spectrograph
**Important** This is a Beta Release Candidate version of this pipeline,
although we have thoroughly tested it there might (and will) still be bugs.
Please let us know if you find any, writting to storres [at] ctio noao edu.

The Goodman High Throughput Spectrograph is an imaging spectrograph,
 if you wish to know more about the instrument please check the 
 [SOAR website](http://www.ctio.noao.edu/soar/content/goodman-high-throughput-spectrograph)
 
To see full documentation please go to the GitHub hosted site for
[Goodman](https://soar-telescope.github.io/goodman/)

## What is contained in this package?

This repository contains tools for processing Goodman's spectrograph data.
It is separated into two main components. The **CCD** part is done by _redccd_
originally developed by [David Sanmartim](https://github.com/dsanmartim)
and currently maintained by our team. It does
the standard _ccd image reduction_, i.e. trim, bias and flat correction.
Currently the ccd reduction pipeline has been integrated into this package so
there is no need to look for it separately. The **Spectroscopic** processing is
done by _redspec_ and includes the following features:


- [x] Identify targets in images
- [x] Trace the target
- [x] Extract the target with background subtraction
- [x] Find the wavelength Solution, interactive and automatically
- [x] Linearize data (resamples)
- [x] Write wavelength solution to a FITS header
- [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different
configurations are provided but is very easy
to add your own.

## How to install it?
### Specific platform instructions
_for installing on Ubuntu 16.04 see 
[this wiki](https://github.com/simontorres/goodman/wiki/Ubuntu-16.04-Installation-Experience)_

_for installing on Centos 7 see 
[this wiki](https://github.com/simontorres/goodman/wiki/Centos-7-Installation)_

Below you will find the instrucctions anyways, please refer to the wiki for the
specifics.

### Get the code

#### Install Git
This step will depend on your platform, below are two examples:

```shell
# Centos
sudo yum install git

# Ubuntu
sudo apt-get install git
```
#### Clone from GitHub
To clone from GitHub copy and paste the following line in a terminal.

```shell
git clone https://github.com/soar-telescope/goodman.git
```
 
You are not ready to install the pipeline yet.

#### Requirements

This software was developed on Python 2.7, use the `requirements.txt` file to 
install all dependencies.

```shell
sudo -H pip2.7 install -r requirements.txt
```

Remember that if you do not have super user privileges, you can still install
all the requirements by adding the `--users` flag. There is also an option of
[installing it using different virtual environments](http://docs.python-guide.org/en/latest/dev/virtualenvs/) 
or even using [Anaconda](https://www.continuum.io/downloads).


You may have some issues when using MatPlotLib on `virtualenvs` and on Mac OSX.
If so, you can try to follow the instructions on
[this site](http://matplotlib.org/faq/osx_framework.html#osxframework-faq) and,
then add the following line on your `.bashrc` or `.bash_profile` file.

```shell
alias python='pythonw' 
```

### Install DCR (Cosmic Ray Rejection)
This pipeline uses [DCR](http://users.camk.edu.pl/pych/DCR/) developed by 
[Wojtek Pych](mailto:pych@camk.edu.pl) instead of `ccdproc.cosmicray_lacosmic` 
because we got better results with `DCR`. Unfortunately you will have to compile
it, I have successfully compiled it on Centos 7, Ubuntu 16.04, Linux Mint 18.1, 
Solaris 11 and MacOS Sierra. The pre-compiled versions are distributed with
the package but it is not guaranteed they will work on your running platform.


Follow [this link](http://users.camk.edu.pl/pych/DCR/) and you can follow the 
instructions there in order to get and compile the `dcr` code. 
The same instructions are reproduced here.

Download the `dcr.tar` file and untar it.
```shell
tar -xvf dcr.tar
```

Compile it
```shell
make
```

If you don't get any errors you can try it without any arguments and you will
get something like this
```shell
$ ./dcr

        USAGE:  dcr  input_file  cleaned_file  cosmicrays_file

File 'dcr.par' must be present in the working directory.
      ~~~~~~
```

#### Make it available for the system
Now that you have compiled the program you have a file called `dcr` you need to
put it in the `$PATH` variable of your system. I usually use `bash` so if you
use something different follow the example below as a guide.

 1. Create a directory to place the executable
     ```shell
     $ mkdir ~/.bin
     ```
     Note that in the example above the directory .bin will be hidden and the 
     symbol `~` denotes your home directory for instance: `/home/goodman/`
 2. Move `dcr` to your new folder.
     ```shell
     $ mv dcr ~/.bin
     ```
 3. Add the directory to the `PATH` variable. With your favorite text editor, 
 open the file `~/.bashrc`
     ```shell
     $ vim ~/.bashrc
     ```
     At the end add the following line.
     ```text
     export PATH=$PATH:/home/user/.bin
     ```
     If you don't know your home directory do the following
     ```shell
     $ cd
     $ pwd
     ```
     Replace `/home/user/`for whatever the output is in the last command
 
 4. Reload the environment variables. For this you can simply close and reopen
 the terminal or you can do:
     ```shell
     $ source ~/.bashrc
     ```

## How to use it?
 
The pipeline is separated in two sub-pipelines. _redccd_ and _redspec_. 
The `--help` argument will print the argument plus some some description

```shell
$ redccd --help
  usage: redccd [-h] [-c] [--ignore-bias] [--auto-clean] [--saturation <Value>]
                [--raw-path raw_path] [--red-path red_path] [--debug]
                [--log-to-file] [--flat-normalize <Normalization Method>]
                [--flat-norm-order <Order>] [--dcr-par-dir <dcr.par directory>]
                [--keep-cosmic-files]

Goodman CCD Reduction - CCDreductions for Goodman spectroscopic data

optional arguments:
  -h, --help            Show this help message and exit
  -c, --cosmic          Clean cosmic rays from science data.
  --ignore-bias         Ignore bias correction
  --auto-clean          Automatically clean reduced data directory
  --saturation <Value>  Saturation limit. Default to 55.000 ADU (counts)
  --raw-path raw_path   Path to raw data.
  --red-path red_path   Path to reduced data.
  --debug               Show detailed information of the process.
  --log-to-file         Write log to a file.
  --flat-normalize <Normalization Method>
                        Choose a method to normalize the master flat
                        forspectroscoy. Choices are: mean, simple (model) and
                        full (fits model to each line).
  --flat-norm-order <Order>
                        Defines the order of the model to be fitted.
  --dcr-par-dir <dcr.par directory>
                        Directory of default dcr.par file.
  --keep-cosmic-files   After cleaning cosmic rays with dcr, do not remove the
                        input file and the cosmic rays file.

  ```

And for `redspec`:

```shell
$ redspec --help
  usage: redspec [-h] [--data-path <Source Path>]
                 [--proc-path <Destination Path>]
                 [--search-pattern <Search Pattern>]
                 [--output-prefix <Out Prefix>] [--extraction <Extraction Type>]
                 [--reference-files <Reference Dir>] [--interactive] [--debug]
                 [--log-to-file] [--save-plots] [--plot-results]

Extracts goodman spectra and does wavelength calibration.

optional arguments:
  -h, --help            show this help message and exit
  --data-path <Source Path>
                        Path for location of raw data. Default <./>
  --proc-path <Destination Path>
                        Path for destination of processed data. Default <./>
  --search-pattern <Search Pattern>
                        Pattern for matching the goodman's reduced data.
  --output-prefix <Out Prefix>
                        Prefix to add to calibrated spectrum.
  --extraction <Extraction Type>
                        Choose a which extraction to perform. Simple is a sum
                        across the spatial direction after the background has
                        been removed. Optimal is a more advanced method that
                        considers weights and profilefitting.
  --reference-files <Reference Dir>
                        Directory of Reference files location
  --interactive         Interactive wavelength solution.Disbled by default.
  --debug               Debugging Mode
  --log-to-file         Write log to a file
  --save-plots          Save all plots in a directory
  --plot-results        Show wavelength calibrated spectrum at the end.

   ```

You should always run `redccd` first and then `redspec`. There are certain
defaults values 

### redccd Defaults

```shell
    --cosmic              False
    --ignore-bias         False
    --auto-clean          False
    --debug               False
    --log-to-file         False
    --keep-cosmic-files   False
    --saturation          55000
    --raw-path            ./
    --red-path            ./RED/
    --flat-normalize      simple
    --dcr-par-dir         files/
    --flat-norm-order     15
```
### redspec Defaults

```shell
    --data-path         ./
    --proc-path         ./
    --search-pattern    cfzsto
    --extraction        simple
    --reference-files   refdata
    --reference-lamp    (empty string)
    --output-prefix     g
    --interactive       False
    --debug             False
    --log-to-file       False
    --save-plots        False
    --plot-results      False
```
