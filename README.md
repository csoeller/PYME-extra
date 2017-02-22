# README #

Plugins and associated code for PYME.

This is currently used as a testbed for

* the new plugin system
* the recipe based processing of SMLM data
* new ideas for data processing

Once sufficiiently mature functionality can migrate into more permanent PYME code if deemed suitable.

### Installation ###

A quick hack so far:

- add the parent directory of this repo to your PYTHONPATH
- add plugin entries in a suitable config file

In my file (.PYME/plugins/visgui/myplugins.txt) I currently have


```
#!python

PYMEcs.experimental.clusterTrack
PYMEcs.experimental.fiducials
PYMEcs.experimental.fiducialsNew

```


### Contact ###

Christian Soeller