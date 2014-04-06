GPU-RSCode
==========
Accelerate Reed-Solomon coding for Fault-Tolerance in RAID-like system.

Online presentation: http://yszheda.github.io/GPU-RSCode/

Motivation
==========
In a RAID-like system, storage is distributed among several devices, and the probability of one of these devices failing becomes significant. Therefore, fault-tolerance must
be taken into account.

Compared to the commonly used replication, erasure codes can reduce the redundancy ratio tremendously. 
Among various types of erasure codes, Reed-Solomon code is one of the popular and it is a kind of optimal erasure codes which has the so-called
Maximum Distance Separable (MDS) property.

The reason that discourages using Reed-Solomon coding to replace replication is its high computational complexity: Galois Field arithmetic is complex,
and matrix operations consume a lot of time. Therefore, GPU acceleratation is taken into account to improve its shortcomings.

Installation
==========
```shell
./configure
make
make install
```

Usage
==========
### New Version ###
```shell
[-h]: show usage information
```
Encode: 
```shell
[-k|-K nativeBlockNum] [-n|-N totalBlockNum] [-e|-E fileName]
```
Decode: 
```shell
[-d|-D] [-k|-K nativeBlockNum] [-n|-N totalBlockNum] 
        [-i|-I originalFileName] [-c|-C config] [-o|-O output]
```
For encoding, the -k, -n, and -e options are all necessary.
For decoding, the -d, -i, and -c options are all necessary.
If the -o option is not set, the original file name will be chosen as the output file name by default.

Performance-tuning Options:
```shell
[-p]: set maxmimum blockDimX
[-s]: set stream number
```
### Old Version (tag 2.0) ###
For encode:
```shell
RS <fragment num> <replica num> -e <original file>
```
For decode:
```shell
RS <fragment num> <replica num> -d [<configuration file>] [-o <output file>]
```
