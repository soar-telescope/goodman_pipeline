# Goodman High Throughput Spectrograph
**Important** This is a Beta Release Candidate version of this pipeline,
although we have thoroughly tested it there might (and will) still be bugs.
Please let us know if you find any.
My e-mail is storres [at] ctio noao edu.

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
- [x] Linearize data (resample)
- [x] Write wavelength solution to a FITS header
- [x] Create a new file for the wavelength calibrated 1D spectrum

There is also a library of calibrated spectrum in FITS format. Different
configurations are provided but is very easy
to add your own.

## How to install it?
### Specific platform instructions
_for installing on Ubuntu 16.04 see [this wiki](https://github.com/simontorres/goodman/wiki/Ubuntu-16.04-Installation-Experience)_

_for installing on Centos 7 see [this wiki](https://github.com/simontorres/goodman/wiki/Centos-7-Installation)_


### Get the code
#### Clone from GitHub
To clone from GitHub copy and paste the following line in a terminal.

```shell
git clone https://github.com/soar-telescope/goodman.git
```
 
You might need to install dependencies, please refer to [Ubuntu](https://github.com/simontorres/goodman/wiki/Ubuntu-16.04-Installation-Experience)
or [Centos7](https://github.com/simontorres/goodman/wiki/Centos-7-Installation) Installation guide. Once you have solved the dependencies issues (if any) do:

```shell
sudo python2.7 setup.py install
```


### Requirements

This software was developed on Python 2.7, use the `requirements.txt` file to install all dependencies.

```shell
sudo -H pip2.7 install -r requirements.txt
```

Remember that if you do not have super user privileges, you can still install all the requirements by adding 
the `--users` flag. There is also an option of [installing it using different virtual environments](http://docs.python-guide.org/en/latest/dev/virtualenvs/) 
or even using [Anaconda](https://www.continuum.io/downloads).


You may have some issues when using MatPlotLib on `virtualenvs` and on Mac OSX. If so, you can try to follow 
the instructions on [this site](http://matplotlib.org/faq/osx_framework.html#osxframework-faq) and, then add the 
following line on your `.bashrc` or `.bash_profile` file.

```bash
alias python='pythonw' 
```


## How to use it?
 
The pipeline is separated in two sub-pipelines. _redccd_ and _redspec_. 
The `--help` argument will print the argument plus some some description

```shell
$ redccd --help
  usage: redccd [-h] [-c] [--ignore-bias] [--auto-clean] [--saturation <Value>]
                [--raw-path raw_path] [--red-path red_path] [--debug]
                [--log-to-file] [--flat-normalize <Normalization Method>]
                [--flat-norm-order <Order>]
  
  Goodman CCD Reduction - CCD reductions for Goodman spectroscopic data
  
  optional arguments:
    -h, --help            show this help message and exit
    -c, --cosmic          Clean cosmic rays from science data.
    --ignore-bias         Ignore bias correction
    --auto-clean          Automatically clean reduced data directory
    --saturation <Value>  Saturation limit. Default to 55.000 ADU (counts)
    --raw-path raw_path   Path to raw data.
    --red-path red_path   Path to reduced data.
    --debug               Show detailed information of the process.
    --log-to-file         Write log to a file.
    --flat-normalize <Normalization Method>
                          Chose a method to normalize the flat for spectroscoy.
                          Choices are: mean, simple (model) and full (fits model
                          to each line).
    --flat-norm-order <Order>
                          Defines the order of the model to be fitted.
  
  ```

And for `redspec`:

   ```shell
     $ redspec --help
     usage: redspec [-h] [-p <Source Path>] [-d <Destination Path>]
                    [-s <Search Pattern>] [-m <Processing Mode>]
                    [-r <Reference Lamp>] [-l <Lamp File>] [-o <Out Prefix>]
                    [-R <Reference Dir>] [-i] [--debug] [--log-to-file]
     
     Extracts goodman spectra and does wavelength calibration.
     
     Supported Processing Modes are:
         <0>: (Default) reads lamps taken at the begining or end of the night.
         <1>: one or more lamps around science exposure.
     
     optional arguments:
       -h, --help            show this help message and exit
       -p <Source Path>, --data-path <Source Path>
                             Path for location of raw data. Default <./>
       -d <Destination Path>, --proc-path <Destination Path>
                             Path for destination of processed data. Default <./>
       -s <Search Pattern>, --search-pattern <Search Pattern>
                             Pattern for matching the goodman's reduced data.
       -m <Processing Mode>, --proc-mode <Processing Mode>
                             Defines the mode of matching lamps to science targets.
       -r <Reference Lamp>, --reference-lamp <Reference Lamp>
                             Name of reference lamp file for mode 0. If not
                             present, the first one in the list will be selected
       -l <Lamp File>        Name of an ASCII file describing which science target
                             uses which lamp. default <lamp.txt>
       -o <Out Prefix>, --output-prefix <Out Prefix>
                             Prefix to add to calibrated spectrum.
       -R <Reference Dir>, --reference-files <Reference Dir>
                             Directory of Reference files location
       -i, --non-interactive
                             Interactive wavelength solution. Enabled by default.
       --debug               Debugging Mode
       --log-to-file         Write log to a file
     
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
    --saturation          55000
    --raw-path            ./
    --red-path            ./RED/
    --flat-normalize      simple
    --flat-norm-order     15
```
### redspec Defaults

```shell
    --data-path         ./
    --proc-path         ./
    --search-pattern    cfzsto
    --proc-mode         0
    --reference-lamp    (empty string)
    --lamp-file         lamps.txt
    --output-prefix     g
    --reference-files   refdata
    --interactive       False
    --debug             False
    --log-to-file       False
```
