# Ctsmr (Continuous Time Stochastic Modelling for R)

This is an example of how to convert a GNU R package to Renjin R.

It is basically an automatic transfer of the `GNU R ctsmr package <http://ctsm.info/>`_ to a Renjin R package 
using the Renjin Maven plugin.

It consists of two modules
- ctsmr-package which converts the GNU R package into a Renjin package
- ctsmr-client which shows a simple usage of the converted package 