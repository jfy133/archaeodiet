#! /usr/bin/env bash

query=$id

raw_res=$(echo "$id" | epost -db taxonomy | esummary -db taxonomy)

parsed_res=$(echo "$raw_res" | sed 's/ /\n/g' | grep -e '<Rank>' -e '<ScientificName>' -e '<CommonName>')

rank==$(echo "$raw_res" | sed 's/ /\n/g' | grep -e '<Rank>' 

if 