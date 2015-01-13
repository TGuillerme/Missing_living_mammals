#!/bin/sh
##########################
#Extract the names of a nexus file and sort them as present in the databases
##########################
#SYNTAX :
#sh TEM_mastsim.sh <nexus matrix>
#with
#<nexus matrix> being the input matrix in nexus format to be sorted
##########################
#version: 0.1
extract.names_version="extract.names.sh v.0.1"
#
#----
#guillert(at)tcd.ie - 13/01/2015
##########################
#Requirements:
#-R 3.x
#-extract.names.R
#-R package "ape"
#-"TaxonReference" folder containing: "FritzTree.rs200k.100trees.tre", "pdb_taxa.csv", and "WilsonReederMSW.csv" 
##########################