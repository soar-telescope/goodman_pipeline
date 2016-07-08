#!/usr/bin/csh


ecl <<END
# < "/home/vultur/iraf/login-xterm.cl"
task \$removekey = /home/vultur/development/SOAR/gmremove.cl

removekey
logout
END
