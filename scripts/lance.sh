#!/bin/bash
set -x
ulimit -S -s unlimited
$1
