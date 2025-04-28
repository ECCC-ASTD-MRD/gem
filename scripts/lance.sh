#!/bin/bash
set -x
ulimit -S -s unlimited
exec $1
