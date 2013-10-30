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
For encode:
```shell
RS <fragment num> <replica num> -c <original file>
```
For decode:
```shell
RS <fragment num> <replica num> -d [<configuration file>]
```

