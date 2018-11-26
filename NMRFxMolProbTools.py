from collections import OrderedDict
import matplotlib.pyplot as plt
import subprocess as sp
import re
import os
import numpy as np


class NMRFxMolProbTools(object):

    # This class field contains some Molprobity parameters we are interested in obversing
    check_list = ('clashscore', 'MolProbityScore', 'pct_badangles',
                  'pct_badbonds', 'ramaOutlier', 'rota<1%', 'numPperpOutliers',
                  'numSuiteOutliers', 'numPperp', 'numSuites')

    def __init__(self):
        self.organized_values = OrderedDict()
        self.file_name = "mpola.txt"

    def org_molprobity_output(self, rfile):
        """ Function reads and organizes output of MolProbity's one line analysis CL function placed in a ".txt" file.
        Returns a dictionary whose keys are the MolProbity tests ran and values are a list of those test
        results/parameters. 
	
	Input Params:

	- rfile (string) : name of the file to read with molprobity output information
	
	Output:
	
	- If file is processed successfully, it'll return an ordered dictionary with molprobity parameters analyzed 
	  as keys and numeric calculations for those parameters as values in a list. 
	  Otherwise, it'll return an empty ordered dictionary.
	"""
	
	# assert the input file is of type 'str'
	assert (type(rfile) is str), "Input file name must be of type: 'str'"

	# check if the provided string is indeed a file
        if os.path.isfile(rfile):

	    # opening and reading file
            with open(rfile, 'r') as f_input:

		# reference pattern for the order of parameters within MolProbity file
                ref_pat = re.compile(r'(\#.+File.*)')

		# MolProbity output information pattern
                analysis_info_pat = re.compile(r"([a-z]*[0-9].*pdb\:.*)")

                file_list = f_input.readlines()
		
                for line in file_list:
                    ref_mat = ref_pat.search(line)
                    analysis_info_mat = analysis_info_pat.search(line)
                    if ref_mat:
                        ref_selection = ref_mat.group(1)
                        ref_split = ref_selection.split(':')
                        for elem in ref_split:
                            self.organized_values[elem] = []

                    if analysis_info_mat:
                        info_selection = analysis_info_mat.group(1)
                        info_split = info_selection.split(':')
                        if self.organized_values:
                            k = 0
                            for _, val in self.organized_values.items():
                                val.append(info_split[k])
                                k += 1

                return self.organized_values

    @classmethod
    def run_descriptive_stats(cls, processed_dict):
        """ This function takes a dictionary and parses the dictionary keys to find specific MolProbity parameters.
        If it finds the parameter, then it will run a series of descriptive statistical calculations and prints
        out those stats.

	"""

        assert (type(processed_dict) is OrderedDict), "Input must be an ordered dictionary."

        if processed_dict:
            # MolProbity Parameters we may want to use :
            # 'clashscore', 'MolProbityScore', 'Mol_pct_rank', 'numbadangles', 'pct_badangles',
            # 'numbadbonds', 'pct_badbonds', 'ramaFavored', 'ramaOutlier', 'rota<1%', 'numPperpOutliers',
            # 'numPperp', 'numSuiteOutliers', 'numSuites'...
            # + more parameters highlighted in MolProb's one-line analysis

            summary = {}
            for key, val_lst in processed_dict.items():
                if key in cls.check_list:
                    try:
                        val = map(float, val_lst)
                    except ValueError:
                        continue
                    if val:
                        val_arr = np.array(val)
                        test_mean = np.mean(val_arr)
                        test_var = np.var(val_arr)
                        test_std = np.std(val_arr)
                        test_min = min(val_arr)
                        test_max = max(val_arr)
                        test_mode = max(set(val), key=val.count)
                        summary[key] = {}
                        summary[key]['mean'] = test_mean
                        summary[key]['var'] = test_var
                        summary[key]['std'] = test_std
                        summary[key]['min'] = test_min
                        summary[key]['max'] = test_max
                        summary[key]['mode'] = test_mode

            # fixme: There has to be an alternative way of calculating these values
            if 'numPperp' and 'numPperpOutliers' in summary:
		# calculate mean and standard deviation outlier % of perpendicular distance from C3'
		# phosphate to the line of glycosidic C1'--N1/9 bond to diagnose strange ribose pucker
                num_Pperp = np.array(map(float,processed_dict['numPperp']))
		num_Pperp_Outliers = np.array(map(float,processed_dict['numPperpOutliers']))
		std_pct_Pperp_Outliers = np.std((num_Pperp_Outliers/num_Pperp)*100)
		mean_pct_Pperp_Outliers = (summary['numPperpOutliers']['mean']/summary['numPperp']['mean'])*100
		
                summary['pct_Pperp'] = {'mean': mean_pct_Pperp_Outliers, 'std': std_pct_Pperp_Outliers}

            if 'numSuites' and 'numSuiteOutliers' in summary:
		# caluculate mean and standard deviation outlier % for suite outliers
		num_Suites = np.array(map(float, processed_dict['numSuites']))
		num_Suite_Outliers = np.array(map(float, processed_dict['numSuiteOutliers']))
		std_pct_Suite_Outliers = np.std((num_Suite_Outliers/num_Suites)*100)
                mean_pct_Suite_Outliers = (summary['numSuiteOutliers']['mean']/summary['numSuites']['mean'])*100

                summary['pct_Suites'] = {'mean': mean_pct_Suite_Outliers, 'std': std_pct_Suite_Outliers}

            return summary


    def run_oneline_analysis(self, final_dir, flag="-dorna -noprotein "):
        """ This function should be used to run a one-line MolProbity analysis on the 'final' directory output
        after running a batch of structures using NMRFxStructure.
        
        Input params:
        
	- final_dir (string): directory of pdb files,
        - flags (string): flags for oneline-analysis command -- initialized for rna analysis.

        Pass an empty string for protein analysis (i.e nmpt.run_oneline_analysis('final', flag=''))
	
	Output:

	- Creates a file with molprobity's oneline-analysis output. Name of file is mpola.txt
	- If a directory 'final_pdbs' doesnt exists, it'll create it and move all the final.pdb files into it

	"""

	# assert the input param is of type 'str'
        assert (type(final_dir) is str), "Input must be of type: 'str'"

	# check if provided path is a directory (oneline-analysis requires a directory containing *.pdb files)
        if os.path.isdir(final_dir):
	    
	    # check if the provided directory contains a final_pdbs directory AND if there are *.pdb files in that directory
            if os.path.isdir(final_dir + "/final_pdbs") and os.path.exists(
                    final_dir + "/final_pdbs/final1.pdb"):

                # name of output file
                ofile_name = final_dir + "/" + self.file_name

		# open output file in 'write' state
		wfile = open(ofile_name, 'w')

		# call subprocess to execute oneline-analysis and write to output file
                out = sp.call(
                    "${HOME}/MolProbity/cmdline/oneline-analysis " + flag + final_dir + "/final_pdbs",
                    shell=True, stdout=wfile)

		wfile.close()
                return

	    # if the above is not true, then check if there are *.pdb files in the provided directory
            if os.path.exists(final_dir + "/final1.pdb"):

		# check if 'final_pdbs' directory exists, if not then create one
                if not os.path.exists(final_dir + "/final_pdbs"):
                    os.mkdir(final_dir + "/final_pdbs")

		    # take all the 'final' pdb files and store it in the 'final_pdbs' directory
                    import glob, shutil

                    try:
                        for pdb_file in glob.glob(final_dir + "/final*.pdb"):
                            shutil.move(pdb_file, final_dir + "/final_pdbs")
                    except:
                        raise ValueError("Something went wrong with *.pdb file transfer process into 'final_pdbs'")

		# output file name
                ofile_name = final_dir + "/" + self.file_name

		# open output file
		wfile = open(ofile_name, 'w')

                out = sp.call(
                    "${HOME}/MolProbity/cmdline/oneline-analysis " + flag + final_dir + "/final_pdbs",
                    shell=True, stdout=wfile)

                wfile.close()
		return
	return None




class ShowMolProbAnalysis(NMRFxMolProbTools):

    def __init__(self, rel_path, file_name="mpola.txt"):
        self.rel_path = rel_path
        self.home = os.getenv("HOME")
        self.file_name = file_name

    def _nav_and_grab(self):
        """Function navigates to the directory and retrieves the output file of MolProbity's one-line analysis. """
	# check if relative path is a string
        if isinstance(self.rel_path, str):
	    # check if relative path exists
            if os.path.exists(self.rel_path):
	    	# check if its a file... if not, check if its a directory, if not 
                if os.path.isfile(self.rel_path):
                    # Pass a file path directly
                    pmp = NMRFxMolProbTools()
                    mp_organized_dict = pmp.org_molprobity_output(self.rel_path)
                    return mp_organized_dict
                elif os.path.isdir(self.rel_path) and os.path.exists(self.rel_path+"/"+self.file_name):
                    # Pass a directory containing the file
                    pmp = NMRFxMolProbTools()
                    mp_organized_dict = pmp.org_molprobity_output(self.rel_path+"/"+self.file_name)
                    return mp_organized_dict
                else:
                    if os.path.isdir(self.rel_path):
                        # Pass a directory containing directories; checks for directories of structures
                        ls_dir = os.listdir(self.rel_path)
                        # PDB ID regex
                        pat = re.compile(r'([1-9][a-z0-9]{3})', re.I)
                        summary_per_struc = []
                        for i, elem in enumerate(ls_dir):
                            match = pat.search(elem)
                            if match:
                                selection = match.group(1)
                                sel_path = self.rel_path+"/"+selection
                                if os.path.isdir(sel_path):
                                    # Found directory of structure in list... now look for file
                                    analysis_file_path = ""
                                    for root, _, afile in os.walk(sel_path):
                                        if self.file_name in afile:
                                            analysis_file_path = root+"/"+self.file_name
                                    if analysis_file_path:
                                        pmp = NMRFxMolProbTools()
                                        processed_dict = pmp.org_molprobity_output(analysis_file_path)
                                        summary = pmp.run_descriptive_stats(processed_dict)
                                        summary_per_struc.append((selection, summary))
                        if summary_per_struc:
                            return summary_per_struc
                        else:
                            raise OSError("File '{}' not found. Check directory or specify 'file_name'.".format(self.file_name))
            else:
                raise ValueError("Could not find relative path used as input. Does not exist.")

    def show_one_struc(self):
        mp_orgnized_dict = self._nav_and_grab()
        assert (type(mp_orgnized_dict) is OrderedDict), ""
        if mp_orgnized_dict:
            obj_arr = ['f{}'.format(file_name.split('.')[0][-1]) for file_name in mp_orgnized_dict['#pdbFileName']]
            index = np.arange(len(obj_arr))
            i = 1
            check_list = NMRFxMolProbTools.check_list
            for k in check_list:
                if k == 'numPperp' or k == 'numSuites':
                    continue
                else:
                    try:
                        value = mp_orgnized_dict[k]
                        if value[0]:
                            value = list(map(float, value))
                        else:
                            continue
                    except KeyError:
                        continue
                    ax = plt.subplot(2, 3, i)
                    ax.bar(index, value, align='edge', width=.5)
                    plt.ylabel(k)
                    plt.ylim(ymin=0)
                    plt.xticks(index, obj_arr)
                    i += 1
            plt.show()

    def show_mult(self):
        """ This function receives a list of tuples whose first element is the structure name (i.e. '2koc')
        and the second element is a dictionary of whose keys are the MolProbity analysis parameter names and whose values
        are dictionaries that summarize descriptive statistical values. """

        # Input list_data example :
        # [('1khm', {'MolProbityScore': {'var': 0.013, 'std': 0.112, 'mean': 2.4516},
        # 'clashscore': {'var': 2.724, 'std': 1.651, 'mean': 4.21},
        # 'numbadbonds': {'var': 0.0, 'std': 0.0, 'mean': 0.0},
        # 'ramaOutlier': {'var': 1.44, 'std': 1.2, 'mean': 5.40},
        # 'rota<1%': {'var': 1.247, 'std': 1.117, 'mean': 3.45},
        # 'numbadangles': {'var': 2.81, 'std': 1.676, 'mean': 5.70}}),
        # ('1z2q', {'MolProbityScore': {'var': 0.027, 'std': 0.163, 'mean': 1.98}, ...
        # 'numbadangles': {'var': 1.087, 'std': 1.043, 'mean': 7.75}})]

        list_data = self._nav_and_grab()

        # Grabbing all of the structure names defined
        pid_object = [name[0] for name in list_data]

        # Defining an index used as a bar plot parameter
        index = np.arange(len(pid_object))

        # Dictionary that pairs the keys of input list_data to the formalize title used in plots
        opts = OrderedDict([
            ('clashscore', 'Clashscore'),
            ('MolProbityScore', 'MolProbity Score'),
            ('rota<1%', 'Rotameric Outliers'),
            ('ramaOutlier', 'Ramanchandran Outliers'),
            ('pct_Suites', 'Suite Outliers (%)'),
            ('pct_Pperp', 'Pucker Outliers (%)')
            # ('pct_badbonds', 'Bad Bonds (%)'),
            # ('pct_badangles', 'Bad Angles (%)')
        ])

        # filtering input list_data to retrieve specified parameters
        info = []
        for key in opts:
            if key in list_data[1][1]:
                info.append((opts.get(key), [tup[1][key]['mean'] for tup in list_data],
                             [tup[1][key]['std'] for tup in list_data]))

        # defining number of plots
        num_of_subplots = len(info)

        if num_of_subplots > 1:
            # helper iterator value
            k = 0
            # iterate over axes and generate subplots based on the length of the 'info' list
            for i in range(1, num_of_subplots + 1):
                ax = plt.subplot((num_of_subplots + 1) // 2, 2, i)
                ax.bar(index, info[k][1], align='center', yerr=info[k][2], color='r',
                       error_kw=dict(ecolor='gray', lw=2, capsize=0))
                plt.ylabel(info[k][0])
                plt.ylim(ymin=0)
                plt.xticks(index, pid_object)
                k += 1

            # Writes title on figure
            plt.suptitle('NMRFxStructure: RNA Structures', fontsize='20')

            plt.show()


# Creating instance
#pmp = NMRFxMolProbTools()
#smp = ShowMolProbAnalysis('../../rnastrucs/', file_name='cy_tests_fa.txt')

#smp.show_mult()
