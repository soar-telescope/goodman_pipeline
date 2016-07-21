#!/usr/bin/csh


ecl <<END
# < "/home/vultur/iraf/login-xterm.cl"
task \$removekey = /home/vultur/development/SOAR/goodman/gmremove.cl

removekey
logout
END

