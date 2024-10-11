In order to use the tools you have to start SPM first so that cat12 is also in the path. 

--------------------------------------------------------------
BA-analysis with some predefined codes:

The tool to try BA is cg_analyze_UKB_long_BA_only. Because reading the huge xls-file takes a lot of time you can set the arguments num and headerwith:
xls_file = '../tables/ukb49261_long3046.xlsx';
[num, header] = cg_get_ukb_data(xls_file);

and then call:
cg_analyze_ukb_long_BAonly(sel, num, header)

The argument 'sel' is the choice of prepared codes. 

--------------------------------------------------------------
VBM or BA-anaylsis:

The function cg_analyze_UKB_long can be also called without argument and you have to define the parameters interactively. Please read the help text for these functions to get info about the arguments.
Always call that function from the analysis folder and organize your data and folders in the following way:
  your_data_folder/analysis
  your_data_folder/tables
  your_data_folder/diff_s6mwmwp1r
  your_data_folder/s6mwmwp1r
  your_data_folder/s6mwmwp2r

You can use the wrapper for cg_analyze_ukb_long with some examples that might help in defining your own analysis. This should be the prefered method, because it is easier to define and log your settings in that way.