#!/bin/sh
#\
exec expect -f "$0" ${1+"$@"}
#-------------------------------------------------------------------------------

# #!/bin/sh
# # \
# exec tclsh "$0" ${1+"$@"}
#-------------------------------------------------------------------------------

#package require Expect

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# QNX: почему-то выходит с ошибкой?
#
#spawn a~.tl TEST_menu @
#ldd:FATAL: Unresolved symbol "openpty" called from libexpect5.45.so
#
# configure:
# checking for openpty... yes
# 
# openpty:
# in libc.a, but not in libc.so
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# in Ubuntu:
# configure:
# checking for openpty... no
# checking for openpty in -lutil... yes
# 
 
#-------------------------------------------------------------------------------
proc expect_press_enter_and_send {spawn_id cmd} {

  expect {
    " ENTER: " {
      send "$cmd\r"
      #puts ""
    }
    timeout {puts "No response from"}
  }

  #  "and PRESS ENTER: " {}

}
#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# puts ""
# puts "CHECK TEST_menu ... "
# puts ""

# #set timeout 2

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # залазим поуправлять программкой с меню (через EXPECT)

# spawn  a~.tl TEST_menu @

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# #expect "and PRESS ENTER: " {puts "\nURRRRAAA!"}
# #expect "FIGNY"

# #log_user 0 ;# останавливает вывод логаута из вызванного процесса !!
# #
# # вообще такой дубовый способ анализа вывода  LOGOUT не очень хорошо,
# # надо бы получать инфу из переменных ... !!

# expect_press_enter_and_send $spawn_id "s"
# expect_press_enter_and_send $spawn_id "3"

# expect_press_enter_and_send $spawn_id "q"
# expect_press_enter_and_send $spawn_id "q"
# # expect_press_enter_and_send $spawn_id "q"

# puts ""
# #puts "~~~~~~~~~~~~~~~~~~~~"

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
proc run_test {ord dx outsteps outfname outnumbs} {

  puts ""
  puts "CHECK TEST ssystem  ... "
  puts ""
  
  spawn ssystem

#  spawn ssystem  2> /dev/null


  expect {
    "Do you want to restore a previously saved integrator state?" {
      send "n\r"
    }
    timeout {puts "No response from 1"}
  }

  expect {
    "Adams order ?" {
      send "$ord\r"
    }
    timeout {puts "No response from 2"}
  }

  expect {
    "step size, days ?" {
      send "$dx\r"
    }
    timeout {puts "No response from 3"}
  }

  expect {
    "Interval between output samples, in steps ?" {
      send "$outsteps\r"
    }
    timeout {puts "No response from 4"}
  }

  expect {
    "Output file name ?" {
      send "$outfname\r"
    }
    timeout {puts "No response from 5"}
  }

  expect {
    "Number of output samples ?" {
      send "$outnumbs\r"
    }
    timeout {puts "No response from 6"}
  }


  # а куда идет дальше вся печать???
  # aaaa, вот ка надо, передать обратно управление:

  interact

}
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# Allocated 2259 bytes to awork
# Allocated 19584 bytes to rwork
# NEQ 102, NEQC 72, 6*NTOTAL 102, IAROIDS 12, NBODY 11, NTOTAL 17
# Reading in physical parameters...
# Do you want to restore a previously saved integrator state?n
# Adams order ? 12
# step size, days ? 0.12500000
# initializing...
# initialized.
# Interval between output samples, in steps ? 320
# Output file name ? test.out
# Number of output samples ? 10

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# set ord  12              ;# Adams order
# set dx   0.125           ;# Step size, days
# set outsteps  320        ;# Interval between output samples
# set outfname "TEST.out"  ;# Output file name
# set outnumbs 10          ;# Number of output samples


if {$argc == 0} {
  puts ""
  puts "Usage: $argv0  <ord> <dx> <outsteps> <outfname> <outnumbs>"
  puts ""
  puts "Example: ssystem.ex 12 0.125 320 TEST.out 10"
  puts ""
  exit
}

set ord       [lindex $argv 0]              
set dx        [lindex $argv 1]     
set outsteps  [lindex $argv 2]       
set outfname  [lindex $argv 3] 
set outnumbs  [lindex $argv 4]         



run_test  $ord $dx $outsteps $outfname $outnumbs

#-------------------------------------------------------------------------------
#
#-------------------------------------------------------------------------------
