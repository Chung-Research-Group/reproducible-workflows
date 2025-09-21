#!/bin/bash

while [ ! -f ./USPEX_IS_DONE ] ; do
   date >> Li4Si1S4-log
   USPEX -r >> Li4Si1S4-log
   sleep 300
done

