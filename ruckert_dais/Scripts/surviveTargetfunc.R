#############################################################
## -file = surviveTargetfunc.R
## -Input: The target range interested in [window] which includes:
##         The lower bound of the range [lbw] &
##         The upper bound of the range [upw]
##
##         The vector of data [sv]
##
## -Output: Which runs in a vector fall within a certain range [Target]
#############################################################
## THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
## NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
## BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
## APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
## AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
##
## -Author: Kelsey Ruckert (klr324@psu.edu)
## -June 10 2015
###############################################################

surviveTarget = function(window, surtarget){
  lbw = window[1]
  upw = window[2]
  sv = surtarget
  
  which((sv >= lbw & sv <= upw))
}
