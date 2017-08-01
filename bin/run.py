import sys
import argparse
from moleculeoperations import *
from argparse import RawTextHelpFormatter

# pymol always wants to get the last word when argparse uses system exit, this stops it
sys.tracebacklimit = 0

parser = argparse.ArgumentParser(description=\
"""-------------------------------------------------------------------------------

=====================
SAAV-GEN-PYMOL-IMAGES
=====================

DESCRIPTION
===========

    This program produces visualizations of single amino
acid variants (SAAVs) mapped onto protein structures, where SAAVs are
represented as small spheres. From a high level, this program takes as its
inputs.

    1) a directory of protein structures for the genomic contexts of interest
       (see --input-dir), i.e. the reference genes that you mapped your
       metagenomic data to. A good question is, how do I get the structure for
       these genes and what if these structures haven't been experimentally
       determined?  The most straightforward way is to use template-based
       structure prediction modelling to "fit" your gene sequence into the
       structure of homologous sequences with known structure. We recommend
       RaptorX Structure Prediction, a free and easy to use service
       (http://raptorx.uchicago.edu/). The accuracy will depend on how similar
       the homologous sequences are to your gene sequence. The RaptorX report
       provides accuracy scores and it is your responsibility to gauge the
       confidence in your predicted structure.
    2) a table of SAAVs identified from metagenomic data (see --saav-table).
       Currently, the only way to identify amino acid variants via metagenomic data is
       with anvio's `anvi-gen-variability-profile`.
    3) a configuration file specifying how the SAAVs are represented (see
       --pymol-config). In this file you specify perspectives and each
       perspective yields a set of images with different settings.

    The program outputs a folder directory (--output-dir) containing the
resultant images and in addition, all of the PyMOL files used to create the
images so they can be explored interactively using the PyMOL interface. After
running this program, a natural next step would be to run
`anvi-saavs-and-protein-structures-summary` to create an HTML output for easy
viewing.

EXAMPLE USAGE
=============

$ pymol -r saav-gen-pymol-images -- --saav-table your_save_table.txt \\
                                    --input-dir your_input_dir \\
                                    --output-dir your_output_dir \\
                                    --pymol-config your_config_file.txt

-------------------------------------------------------------------------------
""", formatter_class=RawTextHelpFormatter)

#==============================================================================

parser.add_argument("-c", "--pymol-config", type=str, required=True, help=\
"""This is the most important parameter of the program. This configuration
file (INI format) provides PyMOL all of the information it needs to
intelligently create those sweet, sweet images you're after. You are able to
control the color, transparency, and sphere_size of each SAAV according to
variables in the SAAV table. Let's look at an example pymol-config file:

    [myperspective1]
    color_var = competing_aas
    color_scheme = competing_aas_cmap
    alpha_var = coverage
    
    [myperspective2]
    color_var = BLOSUM90_weighted
    color_scheme = darkred_to_darkblue
    merged_alpha_var = departure_from_reference
    merged_sphere_size_var = prevalence

In this pymol-config there are two perspectives, therefore two sets of images
will be produced. Perspective names should be descriptive of the options that
define the perspective. The user has cunningly named theirs `myperspective1`
and `myperspective2`.

In `myperspective1`, `color_var = competing_aas` indicates that SAAVs will be
colored according to the `competing_aas` column of the SAAV table (e.g. all
AspGlu's share the same color). `color_scheme = competing_aas_cmap` indicates
the color for each unique value of `competing_aas` is defined according to the
colormap `competing_aas_cmap`, which has been specifically designed for
`color_var = competing_aas` (e.g. AspGlu might map to red and IleVal might map
to blue). For a list of all available color maps, raise an error by providing a
color_scheme that definitely doesn't exist, like `thisshouldntexist`.
`alpha_var = coverage` indicates SAAV tranparency is proportional to the
`coverage` column of the SAAV table. Since the user does not specify a
`sphere_size_var` option, a static sphere_size is chosen by default. If the
user doesn't like the default sphere_size, they can give their own via
`sphere_size_static = x`, where x=2.5, for example.

In `myperspective2`, we see new options prepended with `merged_`. Just how
non-merged options define the images produced for each sample, merged options
define the merged images produced for each group (see --sample-groups for a
definition of groups and merged images). NOTE: Currently only alpha- and
radii-associated options have merged counterparts. With this in mind,
`merged_alpha_var = departure_from_reference` indicates SAAV transparency is
proportional to `departure_from_reference` for all merged images (since
`alpha_var` was not provided, all non-merged images are given a default static
alpha value).  `merged_sphere_size_var = prevalence` is an exceptional case
because `prevalence` is not a column in the SAAV table--but instead of raising
an error, `prevalence` is used as a keyword. So what's prevalence?  It's
defined as the number of samples in a group that have a SAAV at a given codon
position divided by the number of samples in the group.  For example, consider
an AA sequence for 3 samples residing in a group, where "=" means no SAAV and
"o" means SAAV.

    sample1 =o======o=
    sample2 =o======o=
    sample3 =====o==o=
             |   |  |
   position  2   6  9

The prevalences for SAAVs at positions 2, 6, and 9 are respectively, 0.67,
0.33, and 1.00, so in this case SAAV sphere_sizes for merged images are scaled
in proportion to these values. Last thing about `prevalence`: is only a valid
input for `merged_sphere_size_var` or `merged_alpha_var`. The very last thing
to say is that since `merged_color_var` was not provided, it is defaulted to
the value of `color_var`. Given this, we can look back at `myperspective1` and
realize that all `merged_` options defaulted to the non-merged options.
Cool!\n\n""")

#==============================================================================

parser.add_argument("-i", "--input-dir", type=str, required=True, help=\
"""A directory containing .pdb files (.pdb is standard extension for molecular
structures) for the genes of interest. Although the .pdb files may be obtained
through any means imaginable to the user, the directory structure must take the
format of the RaptorX Structure Predicton (http://raptorx.uchicago.edu/)
output, which we recommend using for its simplicity and accuracy when exact
structures are not known. The directory structure should be as follows:

    path/to/input-dir
        |
        |__ <gene_id1>.all_in_one
        |       |
        |       |__ <arbitrary_name1>.pdb
        |       |
        |       |__ ...
        |
        |__ <gene_id2>.all_in_one
        |       |
        |       |__ <arbitrary_name2>.pdb
        |       |
        |       |__ ...
        |
       ...

<gene_id#> should match the gene identifier `corresponding_gene_call` in saav-table.
<arbitrary_name#> can be anything, only the .pdb extension is important. There must
only be one .pdb file in each <gene_id#>.all_in_one folder.\n\n""")

#==============================================================================

parser.add_argument("-t", "--saav-table", type=str, required=True, help=\
"""This is the SAAV table housing all SAAVs and in what genes, and samples they
occur. It is the output of `anvi-gen-variability-profile` with `--engine AA`
envoked.\n\n""")

#==============================================================================

parser.add_argument("-o", "--output-dir", type=str, required=True, help=\
"""The directory image and PyMOL files are output to.\n\n""")

#==============================================================================

parser.add_argument("-g", "--gene-list", type=str, required=False, help=\
"""This is a file containing the genes of interest, identified by values under
the `corresponding_gene_call` of the SAAV table. Here is an example gene-list
file:

    gene_id
    1248
    2342
    2452
    2698

All genes specified must have .pdb files (structure files) in input-dir.\n\n""")

#==============================================================================

parser.add_argument("-s", "--sample-groups", type=str, required=False, help=\
"""Tab-delimited file to specify what samples you are interested in and, if you
want, how these samples group together for the creation of merged images. (If
this file isn't provided, all samples in saav-table will be used and no merged
images will be created). Let's learn through example:

    sample_id   patient ate_breakfast   favorite_color
    170622_1    1   no  blue
    170823_1    1   yes blue
    170622_2    2   no  blue
    170823_2    2   no  blue
    170823_3    3   no  red
    171012_3    3   yes green

The first column must have the header `sample_id` and that's all that's
required for this file to be valid. However, since additional columns were
provided, merged images (protein structures with mapped SAAVs coming from
multiple metagenomes) will be created. For example, merged images for
`ate_breakfast=yes` and `ate_breakfast=no` will be created. (Who else noticed
patient 3 changed his favorite color from red to green? What a dick).\n\n""")

#==============================================================================

parser.add_argument("--ray", type=int, required=False, default=0, help=\
"""Whether or not ray tracing should be performed. Makes images look good but
takes longer""")

#==============================================================================

parser.add_argument("--res", type=int, required=False, default=600, help=\
"""Specifies the width (in pixels) of the output images.""")

#============================================================================== 

args = parser.parse_args()
sys.tracebacklimit = 10

MoleculeOperations(args)

