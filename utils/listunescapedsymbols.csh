#!/bin/sh

nm *.o | grep " T " | egrep -v '(parmetis|METIS|__)' | awk '{printf("#define %s libparmetis__%s\n",$3, $3)}'

