#!/bin/bash

format_clustal() {
    local clustalw_file=$1
    tail -n +2 $clustalw_file | grep -v "^[ *:.]*$"
}



format_clustal test.msa > converted.msa



