"""
TO DO:

DONE 1. scale in proportion to volume, not radius
DONE 2. order color_legend descending to ascending (or alphabetical order)
DONE 3. change raptorXproperty names
DONE 4. add RaptorXProperty variables if found in config but not saav_table
5. config input validity
6. active site analysis
7. sidechain
8. protein_color
9. make color_scheme called `surprisinglythisexists`
DONE 10. add argparse
11. make anvi-append-external-columns-to-SAAV-table binary
"""

import os
import sys
import math
import glob
import time
import pymol
import shutil
import random
import inspect
import functools
import ConfigParser

import numpy as np
import pandas as pd
import pickle as pkl

from pymol import cmd
from colour import Color as colour


__author__ = "Evan Kiefl"
__copyright__ = "Copyright 2017, The anvio Project"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = "1.0"
__maintainer__ = "Evan Kiefl"
__email__ = "kiefl.evan@gmail.com"


class MoleculeOperations():
    
    def __init__(self,  args):
        
    #   define inputs as attributes
        self.args = args
        #A = lambda x: args[x] if x in args else None
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
        self.output_dir = A("output_dir")
        self.input_dir = A("input_dir")
        self.color_vars = A("color_vars")
        self.config_fname = A("pymol_config")
        self.sample_groups_fname = A("sample_groups")
        self.saav_table_fname = A("saav_table")
        self.no_images = A("no_images")
        self.ray = A("ray")
        self.res = A("res")
        self.allow_quince = A("allow_quince")

    #   get dictionary of all columns appended by external classes (for now its just RaptorX)
        self.AddClasses = self.get_dict_of_additional_column_names()

    #   initialize a VariantsTable object
        self.table = VariantsTable(args, rid_quince=True)

    #   make sure the input_dir is in good shape
        self.validate_input_dir(self.input_dir, self.table)

    #   create sample_groups DataFrame
        """ sample_groups.txt very importantly holds any information about how
        the samples may group together. For example, the results from
        clustering algorithms, geographic regions, cohort, or patient groupings """
        self.get_sample_groups()

    #   defines dictionary for how .gif's and .png's are named, based on groupings
        self.get_name_save_dictionary()

    #   load pymol_config file as ConfigParse object and validate its logic
        self.config = Config(self.config_fname, self.AddClasses, self.table)

    #   self.config contains a list of append methods
    #   self.list_of_SAAV_table_append_routines that should be called based on
    #   the values provided in the config file. All of the Classes that append
    #   these variables are then called.
        for Class_ID in self.AddClasses.keys():
            if Class_ID in self.config.list_of_SAAV_table_append_routines:
                self.AddClasses[Class_ID][1](self.table, self.args)

    #   create columns if they don't exist in saav-table already
        add_prevalence = False
    #   read as: if "prevalence" is in your INI file
        for section in self.config.config.sections():
            for option in ["merged_sphere_size_var", "merged_alpha_var"]:
                if self.config.config.has_option(section, option):
                    if self.config.config.get(section, option) == "prevalence":
                        add_prevalence = True
        self.table.merge_sample_groups_to_table(self.sample_groups, add_prevalence)

    #   create the output directory if it doesn't already exist
        self.mkdirp(self.output_dir)

    #   save the config dictionary in self.output_dir as dot file
        self.config.save_pkl(self.output_dir)
    #   save the original sample-groups txt file as a dotfile in self.output_dir
        self.sample_groups.to_csv(os.path.join(self.output_dir, ".sample_groups.txt"), sep='\t', index=False)
    #   save the original gene-list txt file as a dotfile in self.output_dir
        f = open(os.path.join(self.output_dir, ".gene_list.txt"), "w")
        f.write("gene_id\n")
        for gene in self.table.genes:
            f.write("{}\n".format(gene))
        f.close()

    #   loop_through_perspectives calls loop_through_genes which calls loop_through_groupings
        self.loop_through_perspectives()


    def get_dict_of_additional_column_names(self):
        """
        Gets a list of all classes that add columns to the SAAV table from external sources.
        {key : value} looks like {Class_name : ([col1, col2, col3, ...], Class)}
        """
        AddClasses = dict([x for x in inspect.getmembers(sys.modules[__name__], inspect.isclass) if "Add" in x[0]])
        for Class_ID in AddClasses.keys():
            AddClasses[Class_ID] = (AddClasses[Class_ID].column_names(), AddClasses[Class_ID])
        return AddClasses

    def get_name_save_dictionary(self):
        """
        This dictionary defines what is prepended to .pse and .png files and is
        grouping-specific. For example, all groupings except sample_id will
        have "merged" prepended so they are easily parsed and organizable.
        """

        self.name_save = {}
        for grouping in self.sample_groups.columns.values:

            if grouping == "sample_id":
                self.name_save["sample_id"] = "sample_"
            else: 
                self.name_save[grouping] = "merged_{}_".format(grouping)


    def loop_through_perspectives(self):
        """
        1) Makes perspective directory and PyMOL and Images subdirectories
        2) Instantiates colorObject if self.color_hierarchy=="global"
        3) Calls self.loop_through_genes()
        """
        for perspective in self.config.config_dict.keys():

        #   create folder for perspective if doesn't exist
            self.perspective = perspective
            self.perspective_dir = os.path.join(self.output_dir, perspective)
            self.mkdirp(self.perspective_dir)

        #   within the perspective, create three subdirectories: PyMOL, Images, Legends
            self.perspective_dir_pymol = os.path.join(self.perspective_dir, "PyMOL")
            self.mkdirp(self.perspective_dir_pymol)

            self.perspective_dir_legends = os.path.join(self.perspective_dir, "Legends")
            self.mkdirp(self.perspective_dir_legends)
            for t in ["color", "sphere_size", "alpha"]:
                self.mkdirp(os.path.join(self.perspective_dir_legends, t))

            if not self.no_images:
                self.perspective_dir_images = os.path.join(self.perspective_dir, "Images")
                self.mkdirp(self.perspective_dir_images)

        #   This is where the "color_hierarchy" attribute comes in.
        #   colorObject, which determines how the SAAVs are colored, is
        #   reinstantiated based on the color_hierarchy attribute. If
        #   color_hierarchy = global, then colorObject is instantiated only
        #   once (all genes use the same coloring scheme, max/min values, etc.)
        #   If color_hierarchy = gene, then colorObject is reinstantiated for
        #   every gene.  Finally, if color_hierarchy = group, then colorObject
        #   is reinstantiated for every group. The last would be useful if you
        #   you're coloring by competing_aas and there are more than, say, 40
        #   AASTs per gene. To avoid colorObject being reinstantiated
        #   inappropriately, self.get_colormap is always called conditioned by
        #   color_hierarchy being either global, gene, or group. 
            self.color_hierarchy = self.config.config_dict[self.perspective]["color_hierarchy"]
            if self.color_hierarchy == "global":
                saav_table_subset = self.get_relevant_saav_table()
                self.colorObject = Color(saav_table_subset, self.config.config_dict, self.perspective)
                self.colorObject.export_legend(os.path.join(self.perspective_dir_legends, "color"), saav_table_subset, "{}_global_color_legend.txt".format(self.perspective))

        #   loop through each gene
            self.loop_through_genes()


    def loop_through_genes(self):
        """
        1. Makes a gene subdirectory within Images and PyMOL
        2. Generates the protein .pse file
        3. Instantiate colorObject if self.color_hierarchy=="gene"
        4. Call self.loop_through_groupings()
        """
    
        for gene in self.table.genes:

            print("\nCurrently on SECTION: {}, GENE: {}".format(self.perspective, gene))
            
            self.gene = gene
        
        #   get access to protein_pdb
            self.protein_pdb_path = self.get_protein_pdb()

        #   make pymol directory for the gene
            self.gene_dir_pymol = os.path.join(self.perspective_dir_pymol, str(self.gene))
            self.mkdirp(self.gene_dir_pymol)

        #   make image directory for the gene
            if not self.no_images:
                self.gene_dir_images = os.path.join(self.perspective_dir_images, str(self.gene))
                self.mkdirp(self.gene_dir_images)

        #   make the protein .pse 
            self.do_all_protein_pse_things()

        #   create a gene color legend (and reinstantiate colorObject if self.colorhierarchy = gene)
            saav_table_subset = self.get_relevant_saav_table(gene=self.gene)
            if self.color_hierarchy == "gene":
                self.colorObject = Color(saav_table_subset, self.config.config_dict, self.perspective)
            self.colorObject.export_legend(os.path.join(self.perspective_dir_legends, "color"), saav_table_subset, "{}_{}_color_legend.txt".format(self.perspective, self.gene))

        #   loop through each grouping
            self.loop_through_groupings()


    def loop_through_groupings(self):
        """
        Some definitions are needed here for readability.

            groupings : a list of all the groupings, e.g. "sample_id", "cohort", "environment", "region", etc.
            grouping  : a single element of groupings, e.g. "cohort"
            members   : a list of all the categories for a grouping. For example, the grouping "cohort"
                        could have the members "cohort1", "cohort2" and "cohort3"
            member    : a single element of members, i.e. "cohort1"
        """
        
    #   get groupings from sample_groups and loop through them
        self.groupings = self.sample_groups.columns.values
        for grouping in self.groupings:
        
            self.grouping = grouping

        #   for knowing whether merged parameters should be used in the config it must be specified whether grouping is "sample_id" or not.
            self.doing_samples = True if self.grouping == "sample_id" else False

        #   get members from the grouping and loop through them
            self.members = self.sample_groups[self.grouping].unique()
            for member in self.members:

                self.member = member

            #   subset the saav_table to include only group and gene
                saav_table_subset = self.get_relevant_saav_table(gene=self.gene, member=self.member)

            #   redefine alpha and sphere_size max/min (in case alpha_var or sphere_size_var are member-specific (like for prevalence))
                self.sphere_size_and_alpha_min_and_max = self.get_min_and_max_for_alpha_and_sphere_size(saav_table_subset)

            #   make the saav .pse
                self.do_all_saav_pse_things(saav_table_subset)

            #   make an image for the group
                self.do_all_image_things()

            #   delete SAAVs from this group before starting next one
                cmd.delete(self.member)


    def get_min_and_max_for_alpha_and_sphere_size(self, saav_table_subset):
        """
        right now--unlike color--alpha and sphere_size don't get their own class
        or hierarchy parameter so they're defined here. It's cool though.
        """
        sphere_size_and_alpha_min_and_max = {}

        group_specific = ["prevalence"]

        if not self.doing_samples:

            if "merged_alpha_var" in self.config.config_dict[self.perspective].keys():
                if self.config.config_dict[self.perspective]["merged_alpha_var"] in group_specific:
                    sphere_size_and_alpha_min_and_max["merged_alpha_m"] = saav_table_subset[self.grouping+"_"+self.config.config_dict[self.perspective]["merged_alpha_var"]].min()
                    sphere_size_and_alpha_min_and_max["merged_alpha_M"] = saav_table_subset[self.grouping+"_"+self.config.config_dict[self.perspective]["merged_alpha_var"]].max()
                else:
                    sphere_size_and_alpha_min_and_max["merged_alpha_m"] = self.table.saav_table[self.config.config_dict[self.perspective]["merged_alpha_var"]].min()
                    sphere_size_and_alpha_min_and_max["merged_alpha_M"] = self.table.saav_table[self.config.config_dict[self.perspective]["merged_alpha_var"]].max()

            if "merged_sphere_size_var" in self.config.config_dict[self.perspective].keys():
                if self.config.config_dict[self.perspective]["merged_sphere_size_var"] in group_specific:
                    sphere_size_and_alpha_min_and_max["merged_sphere_size_m"] = saav_table_subset[self.grouping+"_"+self.config.config_dict[self.perspective]["merged_sphere_size_var"]].min()
                    sphere_size_and_alpha_min_and_max["merged_sphere_size_M"] = saav_table_subset[self.grouping+"_"+self.config.config_dict[self.perspective]["merged_sphere_size_var"]].max()
                    
                else:
                    sphere_size_and_alpha_min_and_max["merged_sphere_size_m"] = self.table.saav_table[self.config.config_dict[self.perspective]["merged_sphere_size_var"]].min()
                    sphere_size_and_alpha_min_and_max["merged_sphere_size_M"] = self.table.saav_table[self.config.config_dict[self.perspective]["merged_sphere_size_var"]].max()

        else:

            if "sphere_size_var" in self.config.config_dict[self.perspective].keys():
                sphere_size_and_alpha_min_and_max["sphere_size_m"] = self.table.saav_table[self.config.config_dict[self.perspective]["sphere_size_var"]].min()
                sphere_size_and_alpha_min_and_max["sphere_size_M"] = self.table.saav_table[self.config.config_dict[self.perspective]["sphere_size_var"]].max()

            if "alpha_var" in self.config.config_dict[self.perspective].keys():
                sphere_size_and_alpha_min_and_max["alpha_m"] = self.table.saav_table[self.config.config_dict[self.perspective]["alpha_var"]].min()
                sphere_size_and_alpha_min_and_max["alpha_M"] = self.table.saav_table[self.config.config_dict[self.perspective]["alpha_var"]].max()

        return sphere_size_and_alpha_min_and_max


    def do_all_image_things(self):
        """
        """
    #   I could have sophisticated routines here for taking multiple images,
    #   collages, any view setting like orientation, etc. this really deserves
    #   its own class. Instead I'll just call this 1-liner that saves the image
        self.create_image_file() 
        

    def create_image_file(self):
        """
        saves a single png image
        """
        cmd.set_view(self.view)
        if self.ray:
        #   very costly!
            cmd.ray(self.res)

        save_path = os.path.join(self.gene_dir_images, "{}{}.pse".format(self.name_save[self.grouping],self.member))
        cmd.png(save_path)


    def join_pses(self, pse_list):
        """
        Merges multiple pse sessions using the settings of the first file.
    
        INPUT
        -----
        pse_list : list
            list of paths of the pse's to merge. The first one in the list 
            is the one the settings are matched to.
        save : str, default None
            If provided, the merged pse is saved with the filepath `save`
        """
        cmd.load(pse_list[0])
        for pse in pse_list[1:]:
            cmd.load(pse,partial=1)


    def do_all_saav_pse_things(self, saav_table_subset):
        """
        This function creates creates a saav .pse file for each group, and then
        saves the file.
        """
    #   holds all the info for each SAAV like color, sphere_size, & alpha
        self.saav_properties = self.fill_saav_properties_table(saav_table_subset)

    #   create save directory for the saav pse, e.g. "path/to/PyMOL/1248/sample_ANE_132_05M.pse"
        self.saav_pse_path = os.path.join(self.perspective_dir_pymol, str(self.gene), 
                                          "{}{}.pse".format(self.name_save[self.grouping], self.member))

    #   create the file
        self.create_saav_pse_file()


    def create_saav_pse_file(self):
    
    #   set any settings specific to the SAAV .pse files here. Once a SAAV .pse
    #   is merged with its corresponding protein .pse, the settings of the
    #   merged .pse inherits the settings defined under create_protein_pse_file
        cmd.bg_color("white")
        cmd.set("fog","off")

    #   if there are no SAAVs,just save an empty file
        if len(self.saav_properties.index) == 0:
            cmd.save(self.saav_pse_path, self.member)

    #   otherwise do the stuff we were planning to
        else:

        #   define the selection of the SAAVs, create their object, then delete selection
            sites = "+".join([str(resi) for resi in self.saav_properties.index])
            cmd.select(self.member+"_sel","resi {}".format(sites))
            cmd.create(self.member, self.member+"_sel")
            cmd.delete(self.member+"_sel")

            pymol.saav_properties = self.saav_properties

            """ IMPORTANT: In the past PyMOL has not let me perform alter on
            both sphere_transparency and sphere scale without a massive memory
            leak. For now, it seems like this leak is not making itself
            apparent so I have all three alter commands present. Might have to
            fix this. """
        #   change the color for each saav according to the saav_colors dict
            cmd.alter(self.member,"color = pymol.saav_properties.loc[int(resi),'color']")
            cmd.alter(self.member,"s.sphere_transparency = pymol.saav_properties.loc[int(resi),'transparency']")
            cmd.alter(self.member,"s.sphere_scale = pymol.saav_properties.loc[int(resi),'radii']")
            cmd.rebuild()

        #   displays the spheres
            cmd.show("spheres","{} and name ca".format(self.member))

        #   save the file
            cmd.save(self.saav_pse_path, self.member)



    def get_statics_and_variables_for_grouping(self):
        """
        There are lots of keys in the config_dict. This function finds the ones
        that PyMOL cares about for a given grouping. For example, if 
        self.grouping == "sample_id", it might return
        ["color_var", "alpha_static", "sphere_size_var"].
        """
        
        config_keys_for_grouping = [x for x in self.config.config_dict[self.perspective].keys() if \
                                      "_static" in x or "_var" in x]
        if self.doing_samples:
            config_keys_for_grouping = [x for x in config_keys_for_grouping if "merged_" not in x]
        else:
            config_keys_for_grouping = [x for x in config_keys_for_grouping if "merged_" in x]
        return config_keys_for_grouping


    def fill_saav_properties_table(self, saav_table_subset):


    #   if number-type, take the mean. if string-type, take most frequent string (used in the following for loop)
        def resolve_ambiguity(dtype):
            if dtype == "numeric":
                return np.mean
            if dtype == "strings":
                return lambda x: x.value_counts().idxmax()

        def calc_alpha(saav_data_alpha):
            """
            Takes the maximum and minimum values of your alpha data and normalizes it
            to within the range defined by alpha_range (or merged_alpha_range).

            alpha = minimum alpha, beta = maximum alpha, m = mimimum of data, 
            M = maximum of data, a = y-intercept, b = slope
            """

            t = lambda x: x if self.doing_samples else "merged_"+x

            if t("alpha_static") in self.config_keys_for_grouping:
                return saav_data_alpha
            else:
            #   parameters for normalization
                alpha = self.config.config_dict[self.perspective][t("alpha_range")][0]
                beta  = self.config.config_dict[self.perspective][t("alpha_range")][1]
                m = self.sphere_size_and_alpha_min_and_max[t("alpha_m")]
                M = self.sphere_size_and_alpha_min_and_max[t("alpha_M")]
            #   if the max is the min, we don't have a range to normalize over. the best
            #   we can do is return the average of alpha and beta
                if m == M:
                    return (beta + alpha)/2
                a = alpha - m * (beta-alpha) / (M-m)
                b = (beta-alpha) / (M-m)
            #   return normalized version of data
                return a + b * saav_data_alpha


        def calc_sphere_size(saav_data_sphere_size):
            """
            Takes the maximum and minimum values of your sphere_size data and normalizes it
            to within the range defined by sphere_size_range (or merged_sphere_size_range).

            PARAMETERS
            ----------
            alpha, beta :
                minimum, maximum sphere_size
            m, M :
                mimimum, maximum of data
            a, b :
                y-intercep, slopet
            scale_factor :
                scale_factor allows the default values to be user friendly. The default sphere_size_static
                value is 1.0, which PyMOL interprets as a sphere radius of 2. In order for this to be the 
                case, the scale factor ~ 33.5.
            """

            scale_factor = 33.5

            t = lambda x: x if self.doing_samples else "merged_"+x

            if t("sphere_size_static") in self.config_keys_for_grouping:
                return saav_data_sphere_size * scale_factor

            else:
            #   parameters for normalization
                alpha = self.config.config_dict[self.perspective][t("sphere_size_range")][0]
                beta  = self.config.config_dict[self.perspective][t("sphere_size_range")][1]
                m = self.sphere_size_and_alpha_min_and_max[t("sphere_size_m")]
                M = self.sphere_size_and_alpha_min_and_max[t("sphere_size_M")]
            #   if the max is the min, we don't have a range to normalize over. the best
            #   we can do is return the average of alpha and beta
                if m == M:
                    return ( (beta + alpha)/2 ) * scale_factor
                a = alpha - m * (beta-alpha) / (M-m)
                b = (beta-alpha) / (M-m)
            #   return normalized version of data
                return (a + b * saav_data_sphere_size) * scale_factor

    #   first things first, add residue column (+1 used to account for the zero-indexing anvio does)
        saav_table_subset["resi"] = saav_table_subset["codon_order_in_gene"]+1

    #   initialize the two aforementioned DataFrames that use the unique elements of "resi" as their indices
        saav_data = pd.DataFrame({}, index=saav_table_subset["resi"].unique())
        saav_properties = pd.DataFrame({}, index=saav_data.index)

    #   determine for this grouping, the static and variable keys in self.config.config_dict
        self.config_keys_for_grouping = self.get_statics_and_variables_for_grouping()

        methods_for_pymol = {"color"       : self.colorObject.create_color_indices_for_group,
                             "alpha"       : calc_alpha,
                             "sphere_size" : calc_sphere_size}

    #   loop through all the variable parameters
        for key in self.config_keys_for_grouping:

        #   pymol_property is either "color", "alpha", or "sphere_size"
            pymol_property = key.replace("_var","").replace("_static","").replace("merged_","")
            property_value = self.config.config_dict[self.perspective][key]
            if property_value == "prevalence":
                property_value = self.grouping+"_prevalence"

        #   two workflows: one if its variable, and one if its static
            if "_var" in key:

            #   determine whether variable is a string or number datatype. Only color can be string-type
                string_or_number = self.colorObject.color_variable_type if "color_var" in key else "numeric"

            #   add column to saav_data
                saav_data[pymol_property] = \
                    saav_table_subset.groupby("resi").agg({property_value : resolve_ambiguity(string_or_number)})

            elif "_static" in key:

            #   add column to saav_data
                saav_data[pymol_property] = property_value

            else:
                raise ValueError("Dude... You messed up big time.")

            saav_properties[pymol_property] = methods_for_pymol[pymol_property](saav_data[pymol_property])

    #   pymol renders according to transparency, not alpha, so it's transformed here
    #   pymol scales by radius, not sphere size (volume), so it's transformed here
        saav_properties["transparency"] = 1 - saav_properties["alpha"]
        saav_properties["radii"] = ( (3.*saav_properties["sphere_size"]) / (4.*np.pi) ) ** (1./3)
        saav_properties.drop(["alpha","sphere_size"], axis=1, inplace=True)

        return saav_properties


    def do_all_protein_pse_things(self):
        """
        This function creates a directory for the gene if it doesn't exist, locates
        the .pdb file, loads in in a pymol session, applies all programmed settings,
        and saves the file.
        """
    #   loads pdb file into pymol, applies default settings, and saves
        self.create_protein_pse_file()
    #   the protein is in the PyMOL state. now I save it
        self.protein_pse_path = os.path.join(self.perspective_dir_pymol, str(self.gene), "00_{}.pse".format(self.gene))
    #   save it in the gene subfolder and then reinitialize
        cmd.save(self.protein_pse_path)


    def create_protein_pse_file(self):

    #   load and create scaffold and surface objects
        cmd.reinitialize()
        cmd.load(self.protein_pdb_path, "scaffold")
        cmd.copy("surface","scaffold")
        cmd.hide() # creates blank slate to work with

    #   All settings related to the protein .pse should be set here
        cmd.bg_color("white")
        cmd.set("fog", "off")
        cmd.color("wheat", "scaffold")
        cmd.set("ray_trace_mode", "0")
        cmd.set("ray_opaque_background", "off")
        cmd.color("gray90", "surface")
        #cmd.show("surface", "surface")
        cmd.show("cartoon", "scaffold")

        cmd.orient()
        self.view = cmd.get_view()


    def get_sample_groups(self):
        """
        This creates a self.sample_groups DataFrame that specifies

            1) The samples to include
            2) Any clusterings of the samples. 

        Clusterings are used to create composite gifs. sample_groups.txt should
        take the following format:

            sample_id    region    main_groups    proteotypes
            ANE_004_05M    ANE    Gr01    Gr01_C
            ANE_150_05M    ANE    Gr01    Gr01_C
            ANE_151_05M    ANE    Gr01    Gr01_C
            ANE_152_05M    ANE    Gr01    Gr01_C
            ANW_141_05M    ANW    Gr02    Gr02_A
            ANW_142_05M    ANW    Gr02    Gr02_C

        (Afterwards, a method from VariantsTable is envoked to merge the info from self.sample_groups to the
        SAAV table.) Groupings are not required but samples are. If no sample_groups.txt is provided, no 
        groupings are assumed and all the present in the SAAV table are used (after displayinga  warning). 
        Additionally, the group sizes are added to sample groups for prevalence calculations.
        """

    #   if no file provided, accept all sample_ids in SAAV table
        if not self.sample_groups_fname:
            print("\nWARNING: No sample-groups file provided. No groupings will be made and "
                  "samples will be assumed to be those present in the SAAV table.")
            self.sample_groups = pd.Series(self.table.saav_table["sample_id"].unique()).to_frame().rename(columns={0:"sample_id"})

    #   if file was provided but doesn't exist, inform user
        try:
            if not os.path.isfile(self.sample_groups_fname):
                raise ValueError("{} is not a file you moronic piece of scum.".format(self.sample_groups_fname))
        except:
            pass

    #   otherwise, import sample-groups as a pandas DataFrame
        else:
            self.sample_groups = pd.read_csv(self.sample_groups_fname, sep='\t', header=0, index_col=False)

    #   all samples in sample_groups must also be present in the SAAV table
        in_table = list(self.table.saav_table["sample_id"].unique())
        in_samples = list(self.sample_groups["sample_id"].unique())
        in_both = [x for x in in_samples if x in in_table]
        if set(in_both) != set(in_samples):
            raise ValueError("You have samples present in your sample-groups that are not "
                             "in your SAAV table. You can't just do that. The following are in "
                             "sample-groups but not in saav_table".format\
                             ([x for x in in_samples if x not in in_both]))


    def get_relevant_saav_table(self, gene=None, member=None):
        
        saav_table_subset = self.table.saav_table

        if gene:
            saav_table_subset = saav_table_subset[saav_table_subset["corresponding_gene_call"]==gene]
        if member:
            saav_table_subset = saav_table_subset[saav_table_subset[self.grouping]==member]
        return saav_table_subset


    def get_uniques_in_group(self, group_name):
        """
        Returns a list of the unique elements for the column in
        self.sample_groups specified by group_name.  
        """
        return list(self.sample_groups[group_name].unique())


    def get_protein_pdb(self):
        """
        Returns the pdb file for a given gene.
        """ 
        protein_pdbs = glob.glob(os.path.join(self.input_dir, "{}.all_in_one".format(self.gene), "*.pdb"))

        if len(protein_pdbs) != 1:
            raise ValueError("Expecting 1 pdb file but found {}".format(len(protein_pdbs)))
        protein_pdb = protein_pdbs[0]
        return protein_pdb


    def mkdirp(self, path):
        """
        This function makes a directory if it doesn't exist, otherwise it
        doesn't do anything.  It is a python wrapper for the "mkdir -p" command
        in bash
        """
        if os.path.exists(path):
            pass
        else:
            os.mkdir(path)


    """
    Below are a couple of variable validation methods that I've made static so
    they can be borrowed by other classes.
    """

    @staticmethod
    def validate_table(table):
        class_name = table.__class__.__name__
        if not class_name == "VariantsTable":
            raise ValueError(("You have passed an object of class '{}' to the variable "
                              "table. It must be an object of class VariantsTable".\
                              format(class_name)))

    @staticmethod
    def validate_input_dir(input_dir, table):
        """
        checks that the directory exists and that all the gene names match the
        gene names in the saav table
        """
    #   check that folder exists
        if not os.path.isdir(input_dir):
            raise ValueError("you're so bad at this.")
 
    #   if there are zero folders matching the gene_id2.all_in_one format, raise hell
        subdir_list = glob.glob(os.path.join(input_dir, "*.all_in_one"))
        if len(subdir_list) == 0:
            raise ValueError("what the fuck man")
 
    #   if genes in the table aren't found in the raptorx folder, no
        raptor_genes = [int(os.path.splitext(os.path.basename(x))[0]) for x in subdir_list]
        in_both = [gene for gene in table.genes if gene in raptor_genes]
        if not table.genes == in_both:
            missing_in_raptorx = [gene for gene in table.genes if gene not in in_both]
            raise ValueError("You have genes in your table that are missing in your "
                             "structure repository. Here are those that are missing: \n{}".format(missing_in_raptorx))

#======================================================================================

class Color():
    
    def __init__(self, saav_table, pymol_config, perspective):

    #   make input parameters class attributes
        self.saav_table = saav_table
        self.config_dict = pymol_config
        self.perspective = perspective

    #   load HTML color codes and names
    #   (from https://github.com/codebrainz/color-names/blob/master/output/colors.csv)
        self.load_color_db()

    #   This runs when the color is asked to vary
        if "color_var" in self.config_dict[perspective]:
        #   get color_variable from config_dict
            self.color_variable = self.config_dict[perspective]["color_var"]
        #   determine datatype of template_data (is it a string or a number?)
            self.color_variable_type = self.find_color_variable_type()
        #   create template data: an array of the unique entries
        #   if numeric, template_data is just the minimum and maximum value
        #   if string-like, it is all unique entries
            self.template_data = self.get_template_data()
        else:
            self.color_variable = None
            self.color_variable_type = None

    #   if color_static was provided instead of color_var, the color is constant for the perspective
        if "color_static" in self.config_dict[perspective]:
            self.color_static = self.config_dict[perspective]["color_static"]
        else:
            self.color_static = None

    #   construct the colormap differently, depending the datatype of the column
        self.create_colors()

    def access_color(self, x):
        """
        Accesses the color code of x. x is "AspGlu", 2043.20, etc. 
        """

        if self.color_variable_type == "strings":
            return list(self.color_dict[x].rgb)

        if self.color_variable_type == "numeric":
            return list(self.numeric_map(x).rgb)

        if self.color_static:
            return list(colour(self.color_static).rgb)


    def load_color_db(self):
        """
        This loads 866 named colors present on
        https://en.wikipedia.org/wiki/List_of_colors:_A%E2%80%93F and compiled
        into csv format here:
        https://github.com/codebrainz/color-names/blob/master/output/colors.csv
        Useful for string-like data with lots of possible options, for e.g.
        competing AAs. The database is assumed to be in the same directory
        as moleculeoperations.py
        """
        fdir = os.path.dirname(inspect.getfile(self.__class__))
        fname = "colors.csv"
        color_db = pd.read_csv(os.path.join(fdir, fname), names = ("name1","name2","hex","r","g","b"))
        self.color_db = dict(zip(list(color_db["name1"].values), list(color_db["hex"])))

    def create_colors(self):
        """
        This function handles the two types, numeric and strings, differently.
        If strings, a color_dict is made that corresponds all the possible
        inputs to a Color object.  For example, {"AspGlu" : red}. If numeric, a
        color_gradient is made that is a list of Color objects. The numerical
        values are then mapped to these colors by the function
        self.numeric_map. In addition, a color dictionary is made for the
        purposes of exporting a legend, which contains only 3 elements, the
        min, max, and middle data points. For example, if your data runs from
        0-100, the dictionary would look like {0:color1, 50:color2, 100:color3}
        """
        
        if self.color_variable_type == "numeric":
            self.base_colors, self.gradations = self.color_scheme_repo(self.color_variable)
            self.color_gradient = self.get_color_gradient()

        if self.color_variable_type == "strings":
            self.color_dict = self.color_scheme_repo(self.color_variable)
            for key, value in self.color_dict.items():
                self.color_dict[key] = colour(value)

        if self.color_static:
            self.color_dict = {"everything" : colour(self.color_static)}


    def color_scheme_repo(self, var):
        """
        This repo stores default coloring schemes for common variables. These are
        `default` because they can be overwritten by the user. If the variable is
        numeric it returns a base_color list for the construction of a multi-color
        gradient, and if the variable of string-like it directly returns the color
        mapping.
        """

        color_schemes_numeric = {
        "red_to_blue"                 : (["#a03129","#fcf5f4","#ffffff","#e8edf9", "#264799"],[50,25,25,50]),
        "darkred_to_darkblue"         : (["#7c0b03","#7c0b03","#fcf5f4","#ffffff","#e8edf9", "#133382"],[34,50,10,10,50]),
        "white_to_bloodred"           : (["#bc0000","#fceaea"],[150]),
        "rel_diff_from_mean_gene_cov" : (["#b6c7ee","#e8edf9","#ffffff","#fcf5f4", "#7c0b03"],[10,10,10,50])
        }

        color_schemes_string = {
        "ss3_cmap"           : {"alpha helix":"#984ea3", "beta strand":"#ff7f00", "loop":"#1b9e77", "unknown":"#e8e9ea"},
        "solvent_acc_cmap"   : {"buried":"#ff00f2", "intermediate":"#17becf", "exposed":"#004de8", "unknown":"#e8e9ea"}
        }

    #   this is somewhat of a special case since there are too many colors. I just randomly pick 
    #   colors from self.color_db. Otherwise, the color_scheme is chosen from the above dictionaries.
        if var == "competing_aas":
            hex_codes = random.sample(self.color_db.values(), len(self.template_data))
            color_dict = dict(zip(self.template_data, hex_codes))
            return color_dict

    #   numerical types
        if self.color_variable_type == "numeric":
            return color_schemes_numeric[self.config_dict[self.perspective]["color_scheme"]]

    #   string types
        if self.color_variable_type == "strings":
            return color_schemes_string[self.config_dict[self.perspective]["color_scheme"]]


    def get_color_gradient(self):
        """
        This function returns a list of all the colors for a multi-color gradient.
        
        INPUTS
        ------
        self.base_colors : list
            This is a list of strings interpretable by the `colour` module as colors.  
            If you want a 3-color gradient from red to white to blue, then self.base_colors
            should be ["red", "white", "blue"], or [#ff0000, #ffffff, #0000ff].
        self.gradations : integer, list
            The approximate number of self.gradations asked for. If an integer is provided, 
            a color list with len(self.gradations) will be returned (It will not be exactly this
            unless the modulus of gradiations%(n-1) == 0). If you want the number of
            self.gradations to vary between the colors in self.base_colors, you can provide self.gradations
            as a list. In this case, the nth element in self.gradations corresponds to the number
            of self.gradations in the color range defined by the nth and (n+1)th colors in 
            self.base_colors. For example, if self.base_colors = ["red","blue","green"] and self.base_colors
            = [10, 5], 10 colors are defined between red and blue and 5 colors between blue and
            green. The sum of the entries are therefore the total number of color self.gradations
            expected.
        """

        n = len(self.base_colors)

    #   if int-like self.gradations is passed, transform into array-like
        if type(self.gradations) == int:
            self.gradations = [self.gradations//(n-1) for _ in range(n-1)]

    #   if array-like self.gradations is passed, ensure its length is n-1
        if len(self.gradations) != n-1:
                raise ValueError("The self.gradations array must be n-1")

        if any([x==0 for x in self.gradations]):
            raise ValueError("One or more of the gradation intervals you defined sucks.")

        color_gradient = []

        for i in range(1, n):

            fro = colour(self.base_colors[i-1])
            to = colour(self.base_colors[i])

            color_gradient_interval = list(fro.range_to(to, self.gradations[i-1]))
            if not i == n - 1: del color_gradient_interval[-1]
            color_gradient.extend(color_gradient_interval)

        return color_gradient


    def numeric_map(self, x):
        """ 
        This map converts a numeric value to Color object in self.color_gradient.
        """
        n = len(self.color_gradient)

        numeric_range = np.linspace(self.template_data[0], self.template_data[1], n)
        nearest_idx = np.abs(x - numeric_range).argmin()
        return self.color_gradient[nearest_idx]


    def get_template_data(self): 
        """
        Takes SAAV table (already assumed to be subsetted according to
        gene and group) and returns the unique values if string-like and
        returns nothing if numeric type
        """
        if self.color_variable_type == "strings":
            template_data = np.sort(self.saav_table[self.color_variable].unique())
            return template_data

        if self.color_variable_type == "numeric":
            template_data = (self.saav_table[self.color_variable].min(),
                             self.saav_table[self.color_variable].max())
            return template_data

    def find_color_variable_type(self):
        """
        Determines the data type of the column (either string or number)
        """
    #   for some reason strings come up as type == object in pandas
        return "strings" if self.saav_table[self.color_variable].dtype==object else "numeric"

    def create_color_indices_for_group(self, saav_data_color):
        """
        """
        num_saavs = len(saav_data_color.index)
           
        """ IMPORTANT: PyMOL won't let me define color names that contain any
        numbers whatsoever so I have to name the color of each SAAV some 
        unique alphabetic string. This sucks and I'm pissed. """

        alphabet = ["A","B","C","D","E","F","G","H","I","J"]
        color_names = [first + second + third for third in alphabet for second in alphabet for first in alphabet]
        color_names = color_names[:num_saavs]

        color_indices = []
        i = 0
        for resi in saav_data_color.index:
            rgb = self.access_color(saav_data_color.loc[resi])
            cmd.set_color(color_names[i], rgb)
            color_index = [x[1] for x in cmd.get_color_indices() if x[0]==color_names[i]][0]
            color_indices.append(color_index)
            i += 1
        
        return color_indices


    def export_legend(self, path, data, name=None):
        """
        Exports a text file called color_legend.txt
        
        INPUTS
        ------
            path : str
                directory legend is saved to.
            name : name of legend
            data : pandas DataFrame/Series
                data written to the legend
        """

        if not name:
            name = "color_legend.txt"

    #   creates text file
        text_legend = open(os.path.join(path, name), "w")
        text_legend.write("value\tR\tG\tB\thex\n")

        if self.color_static:
        #   get the RGBs as list
            RGB = self.access_color("this argument is useless for this situation :(")
        #   get corresponding hex code
            hexcode = '#%02x%02x%02x' % (RGB[0], RGB[1], RGB[2])
            text_legend.write("{}\t{}\t{}\t{}\t{}\n".format("all", RGB[0], RGB[1], RGB[2], hexcode))

        else:
            
            data = data[self.color_variable]

            if self.color_variable_type == "numeric":
            
            #   if the data type was numeric, we gotta make a color_dict from the data
            #   all this is so the markers aren't ugly (like 4.222222, 8.222224521, etc.)
                data_range = self.template_data[1] - self.template_data[0]
                round_to_what = -int(math.floor(math.log10(data_range/10)))
                markers = np.linspace(self.template_data[0], self.template_data[1], 10)
                markers = np.around(markers, round_to_what)

                for i in range(len(markers)):
                #   get the RGBs as list
                    RGB = self.access_color(markers[i])
                #   colour uses RGB values bounded by [0,1]. I convert to [0,255]
                    RGB = [int(255*X) for X in RGB]
                #   get corresponding hex code
                    hexcode = '#%02x%02x%02x' % (RGB[0], RGB[1], RGB[2])
                #   write this and then on to the next--on, on to the next one
                    text_legend.write("{}\t{}\t{}\t{}\t{}\n".format(markers[i], RGB[0], RGB[1], RGB[2], hexcode))

            else:
                for key in [x for x in sorted(self.color_dict.keys()) if x in list(data.values)]:
                #   get the RGBs as list
                    RGB = self.access_color(key)
                #   colour uses RGB values bounded by [0,1]. I convert to [0,255]
                    RGB = [int(255*X) for X in RGB]
                #   get corresponding hex code
                    hexcode = '#%02x%02x%02x' % (RGB[0], RGB[1], RGB[2])
                #   write this and then on to the next--on, on to the next one
                    text_legend.write("{}\t{}\t{}\t{}\t{}\n".format(key, RGB[0], RGB[1], RGB[2], hexcode))

        text_legend.close()

class AddRaptorXProperty():
    """
    This class holds methods for adding columns to the SAAV table that are
    specifically related to the output of RaptorX Structure Prediction:
    http://raptorx.uchicago.edu/. An input directory is required for the
    instantiation of this class, and should be formatted as follows:

        path/to/input_dir
            gene_id1.all_in_one
                ...
            gene_id2.all_in_one
                ...
            gene_id3.all_in_one
                ...
            ...
    Each folder gene_id#.all_in_one is the standard output from each of the
    RaptorX job submissions.

    If you want to add another method that adds columns to the SAAV table, you
    should add it here and it should start with "append_", then add the added column
    names to the list returned by self.column_names()
    """

    def __init__(self, table, args, ss3_confidence=0.5, ss8_confidence=0.5, solvent_acc_confidence=0.5):

        self.args = args
        #A = lambda x: args[x] if x in args else None
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None

    #   making input variables attributes
        self.input_dir = A("input_dir")
        self.ss3_confidence = ss3_confidence
        self.ss8_confidence = ss8_confidence 
        self.solvent_acc_confidence = solvent_acc_confidence 

    #   make sure table is VariantsTable object
        self.table = table
        MoleculeOperations.validate_table(self.table)

    #   make sure the input_dir is in good shape
        MoleculeOperations.validate_input_dir(self.input_dir, self.table)

    #   define a list of all append methods in this class
        self.all_append_methods = self.get_append_methods()
        
        self.methods_list = self.all_append_methods

    #   add all the columns associated with the append methods in method_list
        self.add_columns()

    @staticmethod
    def column_names():
        return ["ss3",
                "ss3_num_alpha_per_gene",
                "ss3_num_beta_per_gene",
                "ss3_num_loop_per_gene",
                "ss3_num_unknown_per_gene",
                "solvent_acc",
                "solvent_acc_num_buried_per_gene",
                "solvent_acc_num_intermediate_per_gene",
                "solvent_acc_num_exposed_per_gene",
                "solvent_acc_num_unknown_per_gene",
                "ss8",
                "ss8_num_L_per_gene"
                "ss8_num_H_per_gene",
                "ss8_num_G_per_gene",
                "ss8_num_I_per_gene",
                "ss8_num_E_per_gene",
                "ss8_num_B_per_gene",
                "ss8_num_T_per_gene",
                "ss8_num_S_per_gene"]


    def add_columns(self):
        """
        This function calls all the append methods in self.methods_list. If any of the
        methods have optional parameters, they should be passed in the line where
        append is defined. For example, append = Append(parameter = 4)
        """
        for method in self.methods_list:
            append_method = getattr(self, method)
            print("\nAppending method {} from class AddRaptorXProperty".format(method))
            append_method()


    def return_new_table(self):
        """
        Once the new columns have been added, this method should be called by the 
        VariantsTable object.
        """
        return self.table.saav_table


    def get_append_methods(self):
        """
        This function defines a list of all the append methods. It searches for all the methods
        in this class that start with "append_"
        """
        append_methods = [func for func in dir(self) \
                          if callable(getattr(self, func)) and func.startswith("append")] 
        return append_methods


    def append_8_class_structure(self):
        """
        The 8-class secondary structure predictions for RaptorX are found in
        self.input_dir/{}.all_in_one/labeling/{}.ss8. This function uses this
        information from every gene in self.table.saav_table to construct 9 new
        columns:

        ss8 : ["H", "G", "I", "E", "B", "T", "S", "L", "U"]
            The most likely secondary structure. "H" = alpha helix, "E" = beta
            sheet, "C" = loop, ..., and "U" = unknown. We classify the SAAV as
            "U" whenever the highest confidence score attributed to the 3
            structures is less than self.ss8_confidence.
        ss8_genewide_X : integer
            The number of amino acids (not SAAVs--you can count this yourself
            using the table) in the protein that are classified as X
        """
        
        columns = ("codon_order_in_gene","AA","ss8","prob_H","prob_G","prob_I","prob_E","prob_B","prob_T","prob_S","prob_L")

        def calc_ss8_data_in_1gene(x):
        #   get gene id of this groupby object, then find path of .ss8 file for that gene
            gene = x["corresponding_gene_call"].values[0]
            ss8_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.ss8"))[0]
        #   load ss8 data for gene as pandas DataFrame
            ss8 = pd.read_csv(ss8_path, skiprows=3, names=columns, delim_whitespace=True)
        #   0 index the raptor results to conform to anvio convention :\
            ss8["codon_order_in_gene"] -= 1
        #   add gene column (used with 'reference' to uniquely map ss8 entries to saav entries)
            ss8["corresponding_gene_call"] = gene
        #   if the highest confidence for the 8 classes < self.prob_confidence, the SAAV is classified as "U" for unknown
            l = [col for col in ss8.columns if "prob_" in col]
            ss8["ss8"] = ss8.apply(lambda row: row["ss8"] if any(row[l] > self.ss8_confidence) else "U", axis = 1)
        #   ss8_genewide_X is the total number of AAs in the gene with secondary structure X
            ss8["ss8_num_H_per_gene"] = len(ss8[ss8["ss8"] == "H"])
            ss8["ss8_num_G_per_gene"] = len(ss8[ss8["ss8"] == "G"])
            ss8["ss8_num_I_per_gene"] = len(ss8[ss8["ss8"] == "I"])
            ss8["ss8_num_E_per_gene"] = len(ss8[ss8["ss8"] == "E"])
            ss8["ss8_num_B_per_gene"] = len(ss8[ss8["ss8"] == "B"])
            ss8["ss8_num_T_per_gene"] = len(ss8[ss8["ss8"] == "T"])
            ss8["ss8_num_S_per_gene"] = len(ss8[ss8["ss8"] == "S"])
            ss8["ss8_num_L_per_gene"] = len(ss8[ss8["ss8"] == "L"])
            ss8["ss8_num_U_per_gene"] = len(ss8[ss8["ss8"] == "U"])
        #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
        #   ss8_X are also dropped to cut down on saav_table_size, but code could be modified to retain them
            ss8 = ss8.drop(["AA"]+l, axis=1)
        #   merge with original dataframe
            return pd.merge(x,ss8)
        
        saav_table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = saav_table_grouped.apply(calc_ss8_data_in_1gene)


    def append_3_class_structure(self):
        """ 
        The 3-class secondary structure predictions for RaptorX are found in
        self.input_dir/{}.all_in_one/labeling/{}.ss3.  This function uses this
        information from every gene in self.table.saav_table.genes to construct
        4 new columns:

        ss3 : ["H", "E", "C", "U"] The most likely secondary structure. "H" =
        alpha helix, "E" = beta sheet, "C" = loop, and "U" = unknown. We
        classify the SAAV as "U" whenever the highest confidence score
        attributed to the 3 structures is less than self.prob_confidence.

        prob_genewide_X : integer The number of amino acids (not SAAVs--you can
        count this yourself using the table) in the protein that are classified
        as X.  
        """

    #   I am not sure about the names yet, so I added this dictionary so they can change easily
        ss3_rename = {"H" : "alpha helix",
                      "E" : "beta strand",
                      "C" : "loop",
                      "U" : "unknown"}
        
        columns = ("codon_order_in_gene","AA","ss3","prob_H","prob_E","prob_C")

        def calc_ss3_data_in_1gene(x):
        #   get gene id of this groupby object, then find path of .ss3 file for that gene
            gene = x["corresponding_gene_call"].values[0]
            ss3_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.ss3"))[0]
        #   load ss3 data for gene as pandas DataFrame
            ss3 = pd.read_csv(ss3_path, skiprows=3, names=columns, delim_whitespace=True)
        #   0 index the raptor results to conform to anvio convention :\
            ss3["codon_order_in_gene"] -= 1
        #   add gene column (used with 'reference' to uniquely map ss3 entries to saav entries)
            ss3["corresponding_gene_call"] = gene
        #   if the highest confidence for the 3 classes < self.ss3_confidence, the SAAV is classified as "U" for unknown
            l = [col for col in ss3.columns if "prob_" in col]
            ss3["ss3"] = ss3.apply(lambda row: row["ss3"] if any(row[l] > self.ss3_confidence) else "U", axis = 1)
        #   prob_genewide_X is the total number of AAs in the gene with secondary structure X
            ss3["ss3_num_alpha_per_gene"]   = len(ss3[ss3["ss3"] == "H"])
            ss3["ss3_num_beta_per_gene"]    = len(ss3[ss3["ss3"] == "E"])
            ss3["ss3_num_loop_per_gene"]    = len(ss3[ss3["ss3"] == "C"])
            ss3["ss3_num_unknown_per_gene"] = len(ss3[ss3["ss3"] == "U"])
        #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
        #   prob_X are also dropped to cut down on table_size, but code could be modified to retain them
            ss3 = ss3.drop(["AA"]+l, axis=1)
            ss3["ss3"] = ss3["ss3"].replace(to_replace = ss3_rename)
        #   merge with original dataframe
            return pd.merge(x,ss3)

        table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = table_grouped.apply(calc_ss3_data_in_1gene)


    def append_solvent_accessibility(self):
        """
        This function incorporates the solvent accesibility predictions found
        in self.input_dir/{}.all_in_one/labeling/{}.acc by adding one column,
        "solvent_acc".

        solvent_acc : ["B"=Buried(pACC=1-10), "M"=Medium(pACC=11-40),
        "E"=Exposed(pACC=41-100), "U"=Unknown] pACC is equal to the relative
        solvent accessibility calculated by DSSP. If the highest confidence for
        the classifications B, M, and E is less than
        self.solvent_acc_confidence, the SAAV is considered U.
        """

        solvent_acc_rename = {"B" : "buried",
                              "M" : "intermediate",
                              "E" : "exposed",
                              "U" : "unknown"}

        columns = ("codon_order_in_gene","AA","solvent_acc","prob_B","prob_M","prob_E")

        def calc_solv_acc_in_1gene(x):
        #   get gene id of this groupby object, then find path of .ss3 file for that gene
            gene = x["corresponding_gene_call"].values[0]
            solvent_acc_path = glob.glob(os.path.join(self.input_dir,"{}.all_in_one".format(gene),"labeling","*.acc"))[0]
        #   load solvent_acc data for gene as pandas DataFrame
            acc = pd.read_csv(solvent_acc_path, skiprows=3, names=columns, delim_whitespace=True)
        #   0 index the raptor results to conform to anvio convention :\
            acc["codon_order_in_gene"] -= 1
        #   add gene column (used to uniquely map acc entries to saav entries)
            acc["corresponding_gene_call"] = gene
        #   if the highest confidence for the 3 classes < self.solvent_acc_confidence, the SAAV is classified as "U" for unknown
            l = [col for col in acc.columns if "prob_" in col]
            acc["solvent_acc"] = acc.apply(lambda row: row["solvent_acc"] if any(row[l] > self.solvent_acc_confidence) else "U", axis = 1)
        #   prob_genewide_X is the total number of AAs in the gene with secondary structure X
            acc["solvent_acc_num_buried_per_gene"] = len(acc[acc["solvent_acc"] == "B"])
            acc["solvent_acc_num_intermediate_per_gene"] = len(acc[acc["solvent_acc"] == "M"])
            acc["solvent_acc_num_exposed_per_gene"] = len(acc[acc["solvent_acc"] == "E"])
            acc["solvent_acc_num_unknown_per_gene"] = len(acc[acc["solvent_acc"] == "U"])
        #   The `AA` column (1-letter AA code) is redundant, we already have `reference` (3-letter AA code)
        #   prob_X are also dropped to cut down on saav_table_size, but code could be modified to retain them
            acc = acc.drop(["AA"]+l, axis=1)
            acc["solvent_acc"] = acc["solvent_acc"].replace(to_replace = solvent_acc_rename)
        #   merge with original dataframe
            return pd.merge(x,acc)

        saav_table_grouped = self.table.saav_table.groupby("corresponding_gene_call")
        self.table.saav_table = saav_table_grouped.apply(calc_solv_acc_in_1gene)

#==================================================================================================

class VariantsTable():
    
    def __init__(self, args, rid_quince=False):
        
        self.args = args

        #A = lambda x: args[x] if x in args else None
        A = lambda x: args.__dict__[x] if x in args.__dict__ else None
    #   converting input to self variables
        self.saav_table_fname = A("saav_table")
        self.genes_list_fname = A("gene_list")
        self.samples_list_fname = A("samples_list")
        self.simplify_sample_id_method = A("simplify_sample_id_method")
        self.output_dir = A("output_dir")
        self.input_dir = A("raptor_repo")
        self.save_file = A("save_file")

    #   load the saav table
        self.load()

    #   is it quince mode?
        self.quince = True if (self.saav_table["departure_from_consensus"]==0).any() else False
        if self.quince and rid_quince:
            print(("\nWARNING: Your table has `departure_from_consensus` values equal to 0, which probably is the result "
                   "of using the --quince-mode flag during your `anvi-gen-variability-profile` call. This means you were "
                   "about to visualize amino acid positions that are not SAAVs! I am going to save you from this peril "
                   "by removing these entries in your SAAV table."))
            self.saav_table = self.saav_table[self.saav_table["departure_from_consensus"] != 0] 


    #   simplify sample id 
        if self.simplify_sample_id_method:
            self.simplify_sample_ids()

    #   get genes and samples lists
        self.genes = self.load_genes_file_as_list()
        self.samples = self.load_samples_file_as_list()

    #   filter by genes and samples
        self.filter_table("sample_id", self.samples)
        self.filter_table("corresponding_gene_call", self.genes)

    #   get array of columns that exist in this table
        self.columns = self.get_columns()


    def merge_sample_groups_to_table(self, merger, prevalence):
        """
        This function simply merges the information from a sample-groups
        DataFrame into the SAAV table. This is accomplished using the
        pd.merge(df1, df2) method in Pandas. In addition, saav prevalences are
        calculated for each group (e.g. if a saav position is observed in 2
        samples in a group with 10 samples, the prevalence for that saav
        position for that group is 20%.

        INPUTS 
        ------ 
        merger : pandas DataFrame 
            a pandas DataFrame with at least one column name that's found 
            in self.saav_table.
        prevalence : Boolean
            If True, the prevalence scores for each of the groups will be calculated. 
            This is optional because it is slow
        """

    #   add group size information to sample-groups
        def append_group_size(x, group):
            x[group+"_size"] = x["sample_id"].nunique()
            return x

        groups = [group for group in merger.columns if group != "sample_id"]
        for group in groups:
            merger = merger.groupby(group).apply(functools.partial(append_group_size, group=group))

    #   make sure at least one column name is shared
        shared_columns = [x for x in list(merger.columns) if x in list(self.saav_table.columns)]
        if len(shared_columns) == 0:
            raise ValueError("The DataFrame you are trying to merge shares no columns with the SAAV table.")

    #   merge the columns
        self.saav_table = pd.merge(self.saav_table, merger)

    #   now add a prevalence column for each group
        if prevalence:

            def append_group_prevalence(x, group):
                x[group+"_prevalence"] = len(x.index) / x[group+"_size"]
                return x

            groups = [group for group in merger.columns if "_size" not in group and group != "sample_id"]

            for group in groups:
                subset = self.saav_table.groupby(["unique_pos_identifier", group])[group+"_size", "entry_id"].apply(functools.partial(append_group_prevalence, group=group))
                self.saav_table = self.saav_table.merge(subset)

    def get_columns(self):
        """
        Defines an array of the columns in self.saav_table
        """
        return self.saav_table.columns.values


    def load_genes_file_as_list(self):
        """
        Gene file looks like

        gene_id
        1248
        1344
        ...

        """
    #   if there is no file name provided, return all genes in SAAV table
        if self.genes_list_fname is None:
            print(("\nWARNING: You did not provide a list of genes. Including all genes "
                   "present in SAAV table."))
            return list(self.saav_table["corresponding_gene_call"].unique())

    #   otherwise return genes_list
        else:
            if not os.path.isfile(self.genes_list_fname):
                raise ValueError("Your genes-list path sucks, sorry, not sorry.")
            genes = [x.strip() for x in open(self.genes_list_fname).readlines()][1:]
            return [int(x) for x in genes]


    def load_samples_file_as_list(self):
        """
        You really don't know?
        """
    #   if there is no file name provided, return all genes in SAAV table (this might not be all genes in your data!)
        if not self.samples_list_fname:
            return list(self.saav_table["sample_id"].unique())
    #   otherwise return samples_list
        else:
            if not os.path.isfile(self.samples_list_fname):
                raise ValueError("Your samples-list path sucks")
        #   assumes a header
            genes = [x.strip() for x in open(self.samples_list_fname).readlines()][1:]
            return [int(x) for x in genes]


    def filter_table(self, column, elements):
        """
        This function filters the SAAV table according to single column criteria,
        where column is the column you are filtering by and elements is a list of 
        accepted elements in that column. All other rows are filtered
        """
        self.saav_table = self.saav_table[self.saav_table[column].isin(elements)]


    def simplify_sample_ids(self):
        """
        If ad hoc sample_id modifications need to be made for whatever reason, you can define a method
        here for such a purpose. For an example, see the case when simplify_sample_id_method == "SAR11",
        for which this method was first defined.
        """
        if self.simplify_sample_id_method == "SAR11":
            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("SAR11_","")
            self.saav_table["sample_id"] = self.saav_table["sample_id"].str.replace("_BOWTIE2","")
            self.saav_table["cohort"] = self.saav_table["sample_id"].str.split("_",expand=True)[0]
        elif self.simplify_sample_id_method == "placeholder":
            pass
        else:
            raise ValueError("fuck you, man")


    def load(self):
        """
        Load the SAAV table output from gen-variability-profile.
        """
        if not self.saav_table_fname:
            raise ValueError("fuck you... you ass")
        self.saav_table = pd.read_csv(self.saav_table_fname, sep='\t', header=0, index_col=False)


    def save(self):
        """
        Saves any operations performed on the table.
        """
    #   sanity checks
        if not self.save_file:
            raise ValueError("Surely you're joking, Mr. Feynman")
        if os.path.isfile(self.save_file):
            raise ValueError("You can't just DO that...")
    #   save the table in same format style as the output of gen-variability-profile
        self.saav_table.to_csv(self.save_file, sep='\t', index=False)

#==================================================================================================


class Config:
    """
    See help for a description of how the pymol-config file should look.  What
    might be confusing and what's probably even more annoying to explain is the
    merged variables.  All of the non-merged variables define the images
    produced for each sample. But the samples can be grouped according to the
    sample_groups.txt file, and a composite/merged image of each of these
    groupings is also made.  Since many of the SAAVs in the merged image will
    overlap, it's convenient to visualize the degree of overlap with either the
    sphere_size, alpha, or I guess color of the spheres. For this reason, you
    can also define all of the parameters for the merged image. If they are not
    specified, all of the non-merged settings will be copied to the merged
    settings. Okay, that wasn't so bad.
    """
    
    def __init__(self, config_fname, AddClasses, table):

        self.AddClasses = AddClasses
        self.table = table

    #   load the config file
        self.config_fname = config_fname
        self.config = self.load_config_file()

    #   instantiate default values
        self.defaults = self.get_default_dictionary()

    #   all possible options for a given perspective are given here
        self.options_list   = ["color_var",
                               "color_scheme",
                               "color_static",
                               "color_hierarchy",
                               "sphere_size_var",
                               "sphere_size_range",
                               "sphere_size_static",
                               "alpha_var",
                               "alpha_range",
                               "alpha_static",
                               "merged_sphere_size_var",
                               "merged_sphere_size_range",
                               "merged_sphere_size_static",
                               "merged_alpha_var",
                               "merged_alpha_range",
                               "merged_alpha_static",
                               "sidechain"]

    #   initialize list of append routines that may be called based on values provided in pymol-config
        self.list_of_SAAV_table_append_routines = []

    #   try and raise an error
        self.attack_config_structure()

        self.convert_configparser_to_dictionary()

        #self.print_config_dict()


    def save_pkl(self, directory):
        with open(os.path.join(directory, '.config.pkl'), 'wb') as handle:
            pkl.dump(self.config_dict, handle, protocol=pkl.HIGHEST_PROTOCOL)

    def attack_config_structure(self):
        """
        ...make them believe, that offensive operations, often times, is the
           surest, if not the only (in some cases) means of defence - George Washington, 1799
        """

    #   check each perspective name is unique
        if len(self.config.sections()) != len(set(self.config.sections())):
            raise ValueError(("At least one of your perspectives has a name occuring more than once. Well, if we're "
                              "getting technical, at least two of your perspectives have a name occuring more than once."))

    #   check the integrity of each perspective
        for perspective in self.config.sections():

        #   don't allow whitespace in perspective names because why the fuck would you put whitespace in your perspective names?
            if perspective.replace(" ","") is not perspective:
                raise ValueError("Why would you put whitespace in your perspective name? That's not rhetorical. Why? Seriously?")

        #   disallow a pymol_variable (color, alpha, sphere_size) as being both variable and static
            for x in [("color_var", "color_static"),
                      ("sphere_size_var", "sphere_size_static"),
                      ("alpha_var", "alpha_static"),
                      ("merged_sphere_size_var", "merged_sphere_size_static"),
                      ("merged_alpha_var", "merged_alpha_static")]:
                if x[0] in self.config.options(perspective) and x[1] in self.config.options(perspective):
                    raise ValueError(("You cannot specify color, sphere_size, or alpha as being both variable and static. "
                                      "You've specified both {} and {} in perspective {}".format(x[0], x[1], perspective)))

        #   you can't specify a scheme or range if no variable option (..._var) is specified
            for x in [("color_scheme", "color_var"),
                      ("sphere_size_range",  "sphere_size_var"),
                      ("alpha_range",  "alpha_var"),
                      ("merged_sphere_size_range",  "merged_sphere_size_var"),
                      ("merged_alpha_range",  "merged_alpha_var")]:
                if (x[0] in self.config.options(perspective) and x[1] not in self.config.options(perspective)):
                    raise ValueError("You can't specify {} without also specifying {}. (You did so in perspective {})".\
                                      format(x[0], x[1], perspective))

        #   you must specify a color_scheme if variable is selected
            for color_pair in [("color_var", "color_scheme"), ("merged_color_var", "merged_color_scheme")]:
                if color_pair[0] in self.config.options(perspective) and color_pair[1] not in self.config.options(perspective):
                    raise ValueError(("Because you supplied a 'color_var' option in perspective '{}', you must also supply a "
                                     "'color_scheme' option. This error could have been risen due to the lack of a 'merged_color_scheme' "
                                     "option or a 'color_scheme' option. You do the math.".format(perspective)))


            """The composition of options in perspective is sound. Now adding defaults where they are needed"""


        #   if neither _var nor _static options were not provided for any of color, alpha, and sphere_size, static default is set
            for x in [("color_var","color_static"),
                      ("sphere_size_var","sphere_size_static"),
                      ("alpha_var","alpha_static")]:
                if (x[0] not in self.config.options(perspective) and x[1] not in self.config.options(perspective)):
                    self.config.set(perspective, x[1], self.defaults[x[1]])

        #   if color_scheme, sphere_size_range, and alpha_range are not provided (but corresponding _vars are, defaults are set)
            for x in [("color_var","color_scheme"),
                      ("sphere_size_var","sphere_size_range"),
                      ("alpha_var","alpha_range")]:
                if (x[0] in self.config.options(perspective) and x[1] not in self.config.options(perspective)):
                    self.config.set(perspective, x[1], self.defaults[x[1]])

        #   Now for the merged variables. If no merged variables were provided, we copy from the nonmerged variables.
            if len([x for x in self.config.options(perspective) if "merged" in x]) == 0:
                for nonmerged in ["sphere_size_var", "alpha_var", "sphere_size_static", "alpha_static", "sphere_size_range", "alpha_range"]:
                    if nonmerged in self.config.options(perspective):
                        self.config.set(perspective, "merged_"+nonmerged, self.config.get(perspective, nonmerged))

        #   Otherwise, we add defaults the same as for non-merged
            else:
                
            #   if neither _var or _static options were not provided for any of color, alpha, and sphere_size, static default is set
                for x in [("merged_sphere_size_var","merged_sphere_size_static"),
                          ("merged_alpha_var","merged_alpha_static")]:
                    if (x[0] not in self.config.options(perspective) and x[1] not in self.config.options(perspective)):
                        self.config.set(perspective, x[1], self.defaults[x[1]])
            #   if color_scheme, sphere_size_range, and alpha_range are not provided (but corresponding _vars are, defaults are set)
                for x in [("merged_sphere_size_var","merged_sphere_size_range"),
                          ("merged_alpha_var","merged_alpha_range")]:
                    if (x[0] in self.config.options(perspective) and x[1] not in self.config.options(perspective)):
                        self.config.set(perspective, x[1], self.defaults[x[1]])

        #   Finally, we add defaults for variables that are independent of other variables (sidechain, protein color, yada, yada, yada) 
            for option in ["sidechain", "color_hierarchy"]:
                if option not in self.config.options(perspective):
                    self.config.set(perspective, option, self.defaults[option])


            """All defaults have been added. Now checking contents of each option"""
            ############## INCOMPLETE. FOR NOW WE TRUST THE USERS OPTION VALUES ARE VALID ##############


        #   consider each option-value pairing in the 
            for option, value in self.config.items(perspective):

            #   raise hell if user provided options outside of self.options
                if option not in self.options_list:
                    raise ValueError("{} in {} isn't a valid option.".format(option, perspective))

            #   if any of (merged_)alpha_var, color_var, and (merged_)color_var contain columns that are considered extra,
            #   add it to the self.list_of_SAAV_table_append_routines
                if option in ["color_var", "alpha_var", "sphere_size_var", "merged_alpha_var", "merged_sphere_size_var"]:
                    for Class_ID in self.AddClasses.keys():
                        if value in self.AddClasses[Class_ID][0] and value not in list(self.table.saav_table.columns):
                            self.list_of_SAAV_table_append_routines.append(Class_ID)
                            print(("\nWARNING: In perspective '{}' you provided for option '{}' the value '{}', which does "
                                   "not exist in your SAAV table. We will attempt to append this column to your table using "
                                   "the '{}' append class, which at best will take some time, and at worst will fail miserably. "
                                   "Fingers crossed. If you want to create a SAAV table that includes this column permanently, "
                                   "you should run the following command: TBD".\
                                   format(perspective, option, value, Class_ID)))


            """ IMPORTANT """
        #   For now, I don't actually allow merged_color_var,
        #   merged_color_range, or merged_color_scheme to be specified.
        #   Unfortunately the program flow just doesn't easily support
        #   this--it's a design flaw that may be changed in the future. To
        #   foreshadow this potential change, I define the implicit config
        #   variables merged_color_var, merged_color_range, and
        #   merged_color_scheme that the user is unable to define but that
        #   exist implicitly and are simply mirrors of their non-merged
        #   counterparts, color_var, color_scheme, and color_static.
            try:
                self.config.set(perspective, "merged_color_var", self.config.get(perspective, "color_var"))
            except:
                pass
            try:
                self.config.set(perspective, "merged_color_scheme", self.config.get(perspective, "color_scheme"))
            except:
                pass
            try:
                self.config.set(perspective, "merged_color_static", self.config.get(perspective, "color_static"))
            except:
                pass



    def get_default_dictionary(self):
        """
        This dictionary holds the defaults for all the parameters. They are
        only added to perspectives if set explicitly.  What I mean is, if you
        don't provide color_static (for example), this default isn't
        necessarily set. Whether or not it is depends on the logic defined in
        self.attack_config_structure(). For exammple, in the case of
        color_static, it is defined only if color_static and color_var are both
        not specified.
        """

        defaults = {"color_static"              : "#842f68",
                    "color_hierarchy"           : "global",
                    "sphere_size_static"        : "1.0",
                    "merged_sphere_size_static" : "1.0",
                    "alpha_static"              : "0.85",
                    "merged_alpha_static"       : "0.85",
                    "color_scheme"              : "darkred_to_darkblue",
                    "sphere_size_range"         : "0.25, 2.5",
                    "merged_sphere_size_range"  : "0.25, 2.5",
                    "alpha_range"               : "0.10, 0.85",
                    "merged_alpha_range"        : "0.10, 0.85",
                    "sidechain"                 : "no"}
        return defaults


    def load_config_file(self):

    #   check that it exists
        if not os.path.isfile(self.config_fname):
            raise("{} isn't even a file".format(self.config_name))

    #   load as ConfigParser object
        config = ConfigParser.ConfigParser()
        config.read(self.config_fname)

        return config


    def convert_configparser_to_dictionary(self):
        """
        Makes a dictionary out of the ConfigParser object. The reason for
        converting to dictionary is so that we can convert all of the option
        values to the correct data type, e.g. float, boolean, etc. This isn't
        supported by ConfigParser.
        """

        def str_to_bool(x):
            return x.lower() in ("yes", "true", "t", "1")

        def str_to_tuple(x):
            return tuple([float(y) for y in x.replace(" ", "").\
                                              replace("(", "").\
                                              replace(")", "").split(",")])

        type_convert = {"color_var"                 : str,
                        "color_static"              : str,
                        "color_hierarchy"           : str,
                        "merged_color_var"          : str,
                        "merged_color_scheme"       : str,
                        "merged_color_static"       : str,
                        "color_scheme"              : str,
                        "sphere_size_var"           : str,
                        "merged_sphere_size_var"    : str,
                        "sphere_size_static"        : float,
                        "merged_sphere_size_static" : float,
                        "sphere_size_range"         : str_to_tuple,
                        "merged_sphere_size_range"  : str_to_tuple,
                        "alpha_var"                 : str,
                        "merged_alpha_var"          : str,
                        "alpha_static"              : float,
                        "merged_alpha_static"       : float,
                        "alpha_range"               : str_to_tuple,
                        "merged_alpha_range"        : str_to_tuple,
                        "sidechain"                 : str_to_bool}

        self.config_dict = {}
        for perspective in self.config.sections():
            self.config_dict[perspective] = {}
            for option in self.config.options(perspective):
                self.config_dict[perspective][option] = type_convert[option](self.config.get(perspective, option))


    def print_config_dict(self):
        """
        Just some print function
        """
        
        for key in self.config_dict.keys():
            print("{}".format(key))
            for keys, value in self.config_dict[key].items():
                print(keys, value)
            print("\n")


