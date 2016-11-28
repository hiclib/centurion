.. Paris documentation master file, created by
   sphinx-quickstart on Mon Mar 31 17:17:03 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

================================================================================
Centurion: joint identification of centromere locations in yeasts using Hi-C
================================================================================


Although centromeres are essential for life and are the subject of extensive
research, centromere locations in yeast genomes are difficult to infer, and in
most species they are still unknown. Recently, the chromatin conformation
assay Hi-C has been re-purposed for diverse applications, including de novo
genome assembly, deconvolution of metagenomic samples, and inference of
centromere locations. We describe a method, Centurion, that jointly infers the
locations of all centromeres in a single yeast genome by exploiting the
centromeres’ tendency to cluster in 3D space. We first demonstrate the
accuracy of Centurion in identifying known centromere locations from high
coverage Hi-C data of budding yeast and a human malaria parasite. We then use
two metagenomic samples with relatively low coverage Hi-C data to infer
centromere locations for each chromosome in 14 different yeast species. For
yeasts with large centromeres (e.g., S. pombe) Centurion predicts the exact
centromere locations. For seven yeasts with point centromeres, Centurion
predicts most of the centromeres at an average of 5~kb distance from their
known locations. Finally, we predict centromere coordinates for six yeast
species that currently lack centromere annotations. These results suggest that
Centurion can be used for centromere identification for a large number of
yeast species, even with a limited amount of Hi-C sequencing.

Download
========

Download Centurion 0.1 `here
<https://github.com/hiclib/centurion/archive/v0.1.0.tar.gz>`_
or `fork the code on github <https://github.com/hiclib/centurion/>`_.

References
==========

`Accurate identification of centromere locations in yeast genomes using Hi-C
<http://nar.oxfordjournals.org/content/early/2015/05/04/nar.gkv424.full>`_

Contacts
========

If you have any questions or suggestions, please email nelle dot varoquaux at
ensmp dot fr, or open a ticket on `Github
<https://github.com/hiclib/centurion/issues>`_
