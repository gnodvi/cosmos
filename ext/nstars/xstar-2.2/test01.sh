#!/bin/sh
# Activate firefox and paste the clipboard contents into the url bar.

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#wid=`xdotool search --title "Mozilla Firefox"`

# xdotool windowactivate $wid
# sleep 0.2

# xdotool key "ctrl+j"
# xdotool key "BackSpace"
# xdotool key "ctrl+v"
# xdotool key "Return"

#xdotool key "ctrl+t"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

wid=`xdotool search --title "XStar"`

xdotool windowactivate $wid
sleep 0.2
#xdotool key "p"
xdotool key "n"

sleep 3

xdotool windowactivate $wid
sleep 0.2
xdotool key "n"

sleep 3

xdotool windowactivate $wid
sleep 0.2
xdotool key "q"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 
# 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~