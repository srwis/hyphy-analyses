RequireVersion ("2.4.0");
LoadFunctionLibrary("libv3/all-terms.bf");
LoadFunctionLibrary("libv3/UtilityFunctions.bf");
LoadFunctionLibrary("libv3/IOFunctions.bf");
LoadFunctionLibrary("libv3/tasks/estimators.bf");
LoadFunctionLibrary("libv3/tasks/alignments.bf");
LoadFunctionLibrary("libv3/tasks/trees.bf");
LoadFunctionLibrary("libv3/models/codon/MG_REV.bf");
LoadFunctionLibrary("libv3/models/DNA/GTR.bf");
LoadFunctionLibrary("libv3/models/model_functions.bf");
LoadFunctionLibrary("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary("libv3/convenience/math.bf");

LoadFunctionLibrary("SelectionAnalyses/modules/io_functions.ibf");
LoadFunctionLibrary("SelectionAnalyses/modules/selection_lib.ibf");

utility.SetEnvVariable ("NORMALIZE_SEQUENCE_NAMES", TRUE);

/*----------Branch-Site Directional Selection Analysis---------*/
/*Usage: Takes an alignment and a tree (with foreground branches tagged with "{FG}"), a nucleotide model
string ("010010" for HKY85, "012345" for REV), a range of sites to test, and a site reference shift,
which allows for alignments that don't start where they should*/

meds.analysis_description = {terms.io.info : "RUNNING MEDS (Models for Episodic Directional Selection). For help please refer to https://github.com/veg/hyphy-analyses .",
                               terms.io.version : "1.1",
                               terms.io.authors : "Ben Murrell, Sergei L Kosakovsky Pond, Sadie Wisotsky",
                               terms.io.contact : "spond@temple.edu",
                               terms.io.requirements : "codon alignment and a phylogenetic tree"
                              };


/*--Analysis Setup--*/

io.DisplayAnalysisBanner (meds.analysis_description);
//fprintf (stdout, "\n[RUNNING MEDS (Models for Episodic Directional Selection). For help please refer to https://github.com/veg/hyphy-analyses]\n");

meds.json    = { terms.json.analysis: meds.analysis_description,
					 terms.json.input: {},
					 terms.json.fits : {},
					 terms.json.timers : {},
					 meds.terms.json.site_logl : {},
					 meds.terms.json.evidence_ratios : {},
					 meds.terms.json.site_reports : {}
					};

KeywordArgument ("code",        "Which genetic code should be used", "Universal");
KeywordArgument ("alignment",   "An codon alignment in one of the formats supported by HyPhy");
KeywordArgument ("tree",    "Please select a tree file for the data:");
KeywordArgument ("branches", "Select Branches", "All");

namespace meds {
    LoadFunctionLibrary ("SelectionAnalyses/modules/shared-load-file.bf");
    load_file ({utility.getGlobalValue("terms.prefix"): "meds"});
}

meds.table_screen_output  = {{"Codon", "Partition", "alpha", "beta", "LRT", "Selection detected?"}};
meds.table_output_options = {terms.table_options.header : TRUE, terms.table_options.minimum_column_width: 16, terms.table_options.align : "center"};


//just setting this to right Now
meds.srv = FALSE;
meds.site_alpha = "Site relative synonymous rate";
meds.site_beta = "Site relative non-synonymous rate (tested branches)";
meds.site_beta_nuisance = "Site relative non-synonymous rate (untested branches)";

// default cutoff for printing to screen
meds.p_value = 0.1;


//unsure if this is needed
meds.branch_sets = {};

utility.ForEachPair (meds.selected_branches[0], "_branch_", "_model_",
"
    utility.EnsureKey (meds.branch_sets, _model_);
    meds.branch_sets[_model_] + _branch_;
");

meds.branch_class_count = utility.Array1D (meds.branch_sets);



/** this will fit the GTR model to the data read previously
    and populate the dictionary of results from the fit as
    namespace.gtr_results (namespace is the argument you pass
    to doGTR)
*/

namespace meds {
    doGTR ("meds");
}

/** this is a standard library call to mark all global parameters of this fit as
fixed for subsequent fits; this is accomplished by adding a "fix-me" (terms.fix)
field to the parameter record  */

estimators.fixSubsetOfEstimates(meds.gtr_results, meds.gtr_results[terms.global]);



/** this will the MG94xREV model with multiple branch classes, each having its own
    dN/dS. Branch lengths are proportional to those obtained with the GTR fit.
    This will create 'namespace'.partitioned_mg_results
*/

namespace meds {
    scaler_prefix = "meds.scaler";
    doPartitionedMG ("meds", FALSE);
}

io.ReportProgressMessageMD ("MEDS", "codon-refit", "Improving branch lengths, nucleotide substitution biases, and global dN/dS ratios under a full codon model");

meds.final_partitioned_mg_results = estimators.FitMGREV (meds.filter_names, meds.trees, meds.codon_data_info [terms.code], {
    terms.run_options.model_type: terms.local,
    terms.run_options.partitioned_omega: meds.selected_branches,
    terms.run_options.retain_lf_object: TRUE
}, meds.partitioned_mg_results);


io.ReportProgressMessageMD("MEDS", "codon-refit", "* " + selection.io.report_fit (meds.final_partitioned_mg_results, 0, meds.codon_data_info[terms.data.sample_size]));
meds.global_dnds = selection.io.extract_global_MLE_re (meds.final_partitioned_mg_results, "^" + terms.parameters.omega_ratio);

utility.ForEach (meds.global_dnds, "_value_", 'io.ReportProgressMessageMD ("MEDS", "codon-refit", "* " + _value_[terms.description] + " = " + Format (_value_[terms.fit.MLE],8,4));');

// define the site-level likelihood function

meds.site.mg_rev = model.generic.DefineModel("models.codon.MG_REV.ModelDescription",
        "meds_mg", {
            "0": parameters.Quote(terms.local),
            "1": meds.codon_data_info[terms.code]
        },
        meds.filter_names,
        None);


meds.site_model_mapping = {"meds_mg" : meds.site.mg_rev};



/* set up the local constraint model */

meds.alpha = model.generic.GetLocalParameter (meds.site.mg_rev, utility.getGlobalValue("terms.parameters.synonymous_rate"));
meds.beta = model.generic.GetLocalParameter (meds.site.mg_rev, utility.getGlobalValue("terms.parameters.nonsynonymous_rate"));
io.CheckAssertion ("None!=meds.alpha && None!=meds.beta", "Could not find expected local synonymous and non-synonymous rate parameters in \`estimators.FitMGREV\`");

selection.io.startTimer (meds.json [terms.json.timers], "MEDS analysis", 2);

//----------------------------------------------------------------------------------------

// I think this function scales the branches by alpha and beta for each site????
function meds.apply_proportional_site_constraint (tree_name, node_name, alpha_parameter, beta_parameter, alpha_factor, beta_factor, branch_length) {

    meds.branch_length = (branch_length[terms.parameters.synonymous_rate])[terms.fit.MLE];

    node_name = tree_name + "." + node_name;

    ExecuteCommands ("
        `node_name`.`alpha_parameter` := (`alpha_factor`) * meds.branch_length__;
        `node_name`.`beta_parameter`  := (`beta_factor`)  * meds.branch_length__;
    ");
}
//----------------------------------------------------------------------------------------
//adds scalers to the global parameters
meds.scalers = {{"meds.alpha_scaler", "meds.beta_scaler_test", "meds.beta_scaler_nuisance"}};

model.generic.AddGlobal (meds.site.mg_rev, "meds.alpha_scaler", meds.site_alpha);
model.generic.AddGlobal (meds.site.mg_rev, "meds.beta_scaler_test", meds.site_beta);
model.generic.AddGlobal (meds.site.mg_rev, "meds.beta_scaler_nuisance", meds.site_beta_nuisance);
parameters.DeclareGlobal (meds.scalers, {});


//this is the stuff I think I need from meds???
lfunction meds.rate.modifier (fromChar, toChar, namespace, model_type, model) {
    baseline = Call (^"meds.site.mg_rev.rate", fromChar,toChar, namespace, model_type, model);
    utility.EnsureKey (baseline, model_type);
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.meds.bias"), utility.getGlobalValue("meds.parameter.bias"));
    selection.io.json_store_key_value_pair (baseline, model_type, utility.getGlobalValue("terms.meds.rate"), utility.getGlobalValue("meds.parameter.rate"));
    baseline [utility.getGlobalValue("terms.model.rate_entry")] = parameters.AppendMultiplicativeTerm (baseline [utility.getGlobalValue("terms.model.rate_entry")], utility.getGlobalValue("meds.parameter.rate"));
    if (toChar == model["meds.residue_bias"]) {
        baseline [utility.getGlobalValue("terms.model.rate_entry")] =
            parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                  "`utility.getGlobalValue("meds.parameter.bias")`/(1-Exp (-`utility.getGlobalValue("meds.parameter.bias")`))");
     } else {
        if (fromChar == model["meds.residue_bias"]) {
            parameters.AppendMultiplicativeTerm ( baseline [utility.getGlobalValue("terms.model.rate_entry")],
                                                  "`utility.getGlobalValue("meds.parameter.bias")`/(Exp (`utility.getGlobalValue("meds.parameter.bias")`-1))");
        }
    }
    return baseline;
}
lfunction meds.biased.model.generator (type, residue) {
    model = Call (^"meds.site.mg_rev", type);
    utility.setGlobalValue("meds.site.mg_rev.rate", model[utility.getGlobalValue ("terms.model.q_ij")]);
    model[utility.getGlobalValue ("terms.model.q_ij")] = "meds.rate.modifier";
    model["meds.residue_bias"] = residue;
    return model;
}

//----------------------------------------------------------------------------------------

// I think this takes a likelihood and data and applies the scalers to each site
// then tests for selection using the alt model where everything is allowed to vary?
// and the null where everything is 1????
lfunction meds.handle_a_site (lf, filter_data, partition_index, pattern_info, model_mapping) {



    GetString (lfInfo, ^lf,-1);
    ExecuteCommands (filter_data);
    __make_filter ((lfInfo["Datafilters"])[0]);

    utility.SetEnvVariable ("USE_LAST_RESULTS", TRUE);

    if (^"meds.srv"){
        ^"meds.alpha_scaler" = 1;
    } else
    {
        ^"meds.alpha_scaler" := 1;
    }
    ^"meds.beta_scaler_test"  = 1;
    ^"meds.beta_scaler_nuisance"  = 1;

    Optimize (results, ^lf);
    console.log( "this");

    alternative = estimators.ExtractMLEs (lf, model_mapping);
    alternative [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];

    ^"meds.alpha_scaler" = (^"meds.alpha_scaler" + 3*^"meds.beta_scaler_test")/4;
    parameters.SetConstraint ("meds.beta_scaler_test","meds.alpha_scaler", "");

    Optimize (results, ^lf);
    console.log( "that");

    Null = estimators.ExtractMLEs (lf, model_mapping);


    Null [utility.getGlobalValue("terms.fit.log_likelihood")] = results[1][0];


    Export (lfs, ^lf);
    fprintf (MESSAGE_LOG, lfs);
    //assert (0);

    return {utility.getGlobalValue("terms.alternative") : alternative, utility.getGlobalValue("terms.Null"): Null};
}

/* echo to screen calls */
meds.report.counts        = {{0,0}};



/* I don't think I need this bit for MEDS
meds.report.positive_site = {{"" + (1+((meds.filter_specification[meds.report.partition])[terms.data.coverage])[meds.report.site]),
                                    meds.report.partition + 1,
                                    Format(meds.report.row[0],10,3),
                                    Format(meds.report.row[1],10,3),
                                    Format(meds.report.row[3],10,3),
                                    "Pos. p = " + Format(meds.report.row[4],6,4)}};

meds.report.negative_site = {{"" + (1+((meds.filter_specification[meds.report.partition])[terms.data.coverage])[meds.report.site]),
                                    meds.report.partition + 1,
                                    Format(meds.report.row[0],10,3),
                                    Format(meds.report.row[1],10,3),
                                    Format(meds.report.row[3],10,3),
                                    "Neg. p = " + Format(meds.report.row[4],6,4)}};
*/

//think this is how everything get spit out into the markdown format???
//don't think I need it either
/*
function meds.report.echo (meds.report.site, meds.report.partition, meds.report.row) {

    meds.print_row = None;
    if (meds.report.row [4] < meds.pvalue) {
        if (meds.report.row[0] < meds.report.row[1]) {
            meds.print_row = meds.report.positive_site;
            meds.report.counts[0] += 1;
        } else {
            meds.print_row = meds.report.negative_site;
            meds.report.counts [1] += 1;
        }
    }

     if (None != meds.print_row) {
            if (!meds.report.header_done) {
                io.ReportProgressMessageMD("meds", "" + meds.report.partition, "For partition " + (meds.report.partition+1) + " these sites are significant at p <=" + meds.pvalue + "\n");
                fprintf (stdout,
                    io.FormatTableRow (meds.table_screen_output,meds.table_output_options));
                meds.report.header_done = TRUE;
                meds.table_output_options[terms.table_options.header] = FALSE;
            }

            fprintf (stdout,
                io.FormatTableRow (meds.print_row,meds.table_output_options));

        }

}
*/


//do need results
lfunction meds.store_results (node, result, arguments) {

    partition_index = arguments [2];
    pattern_info    = arguments [3];

    result_row          = { { 0, // alpha
                          0, // beta
                          0, // alpha==beta
                          0, // LRT
                          1, // p-value,
                          0  // total branch length of tested branches
                      } };


    if (None != result) { // not a constant site

        lrt = math.DoLRT ((result[utility.getGlobalValue("terms.Null")])[utility.getGlobalValue("terms.fit.log_likelihood")],
                          (result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.fit.log_likelihood")],
                          1);
        result_row [0] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"meds.site_alpha");
        result_row [1] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.alternative")], ^"meds.site_beta");
        result_row [2] = estimators.GetGlobalMLE (result[utility.getGlobalValue("terms.Null")], ^"meds.site_beta");
        result_row [3] = lrt [utility.getGlobalValue("terms.LRT")];
        result_row [4] = lrt [utility.getGlobalValue("terms.p_value")];



        sum = 0;
        alternative_lengths = ((result[utility.getGlobalValue("terms.alternative")])[utility.getGlobalValue("terms.branch_length")])[0];

        utility.ForEach (^"meds.case_respecting_node_names", "_node_",
                '_node_class_ = ((^"meds.selected_branches")[`&partition_index`])[_node_];
                 if (_node_class_ == utility.getGlobalValue("terms.tree_attributes.test")) {
                    `&sum` += ((`&alternative_lengths`)[_node_])[utility.getGlobalValue("terms.fit.MLE")];
                 }
            ');

        result_row [5] = sum;
    }


    utility.EnsureKey (^"meds.site_results", partition_index);
/* I don't think I need this
    utility.ForEach (pattern_info[utility.getGlobalValue("terms.data.sites")], "_meds_result_",
        '
            (meds.site_results[`&partition_index`])[_meds_result_] = `&result_row`;
            meds.report.echo (_meds_result_, `&partition_index`, `&result_row`);
        '
    );

*/
    //assert (0);
}
//----------------------------------------------------------------------------------------

meds.site_results = {};

for (meds.partition_index = 0; meds.partition_index < meds.partition_count; meds.partition_index += 1) {
    meds.report.header_done = FALSE;
    meds.table_output_options[terms.table_options.header] = TRUE;
    model.ApplyModelToTree( "meds.site_tree", meds.trees[meds.partition_index], {terms.default : meds.site.mg_rev}, None);

    meds.case_respecting_node_names = trees.branch_names (meds.site_tree, TRUE);
    //console.log(med.case_respecting_node_names);

    meds.site_patterns = alignments.Extract_site_patterns ((meds.filter_specification[meds.partition_index])[utility.getGlobalValue("terms.data.name")]);

    // apply constraints to the site tree
    // alpha = alpha_scaler * branch_length
    // beta  = beta_scaler_test * branch_length or beta_nuisance_test * branch_length

    utility.ForEach (meds.case_respecting_node_names, "_node_",
            '_node_class_ = (meds.selected_branches[meds.partition_index])[_node_];
             if (_node_class_ == terms.tree_attributes.test) {
                _beta_scaler = meds.scalers[1];
             } else {
                _beta_scaler = meds.scalers[2];
             }
             meds.apply_proportional_site_constraint ("meds.site_tree", _node_, meds.alpha, meds.beta, meds.scalers[0], _beta_scaler, (( meds.final_partitioned_mg_results[terms.branch_length])[meds.partition_index])[_node_]);
        ');


    // create the likelihood function for this site
    ExecuteCommands (alignments.serialize_site_filter
                                       ((meds.filter_specification[meds.partition_index])[utility.getGlobalValue("terms.data.name")],
                                       ((meds.site_patterns[0])[utility.getGlobalValue("terms.data.sites")])[0],
                   ));

    __make_filter ("meds.site_filter");

    LikelihoodFunction meds.site_likelihood = (meds.site_filter, meds.site_tree);



    estimators.ApplyExistingEstimates ("meds.site_likelihood", meds.site_model_mapping, meds.final_partitioned_mg_results,
                                        terms.globals_only);


    meds.queue = mpi.CreateQueue ({"LikelihoodFunctions": {{"meds.site_likelihood"}},
                                   "Models" : {{"meds.site.mg_rev"}},
                                   "Headers" : {{"libv3/all-terms.bf"}},
                                   "Variables" : {{"meds.srv"}}
                                 });


    /* run the main loop over all unique site pattern combinations */
    utility.ForEachPair (meds.site_patterns, "_pattern_", "_pattern_info_",
        '
            if (_pattern_info_[terms.data.is_constant]) {
                meds.store_results (-1,None,{"0" : "meds.site_likelihood",
                                                                "1" : None,
                                                                "2" : meds.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : meds.site_model_mapping
                                                                });
            } else {
                mpi.QueueJob (meds.queue, "meds.handle_a_site", {"0" : "meds.site_likelihood",
                                                                "1" : alignments.serialize_site_filter
                                                                   ((meds.filter_specification[meds.partition_index])[terms.data.name],
                                                                   (_pattern_info_[utility.getGlobalValue("terms.data.sites")])[0]),
                                                                "2" : meds.partition_index,
                                                                "3" : _pattern_info_,
                                                                "4" : meds.site_model_mapping
                                                                },
                                                                "meds.store_results");
            }
        '
    );

    mpi.QueueComplete (meds.queue);
    meds.partition_matrix = {Abs (meds.site_results[meds.partition_index]), Rows (meds.table_headers)};

    utility.ForEachPair (meds.site_results[meds.partition_index], "_key_", "_value_",
    '
        for (meds.index = 0; meds.index < Rows (meds.table_headers); meds.index += 1) {
            meds.partition_matrix [0+_key_][meds.index] = _value_[meds.index];
        }
    '
    );

    meds.site_results[meds.partition_index] = meds.partition_matrix;


}

meds.json [terms.json.MLE ] = {terms.json.headers   : meds.table_headers,
                               terms.json.content : meds.site_results };


io.ReportProgressMessageMD ("meds", "results", "** Found _" + meds.report.counts[0] + "_ sites under pervasive positive diversifying and _" + meds.report.counts[1] + "_ sites under negative selection at p <= " + meds.pvalue + "**");


return 0;


/*---------Loading alignment and tree files-------------------------------------*/


//fprintf (stdout, "Loaded ", myData.species, " sequences with ", myData.sites, " sites from ",LAST_FILE_PATH,"\n");


KeywordArgument ("output", "Write", "test.csv");
//SetDialogPrompt ("Specify the output (.csv) file");


siteShift 		= -1;  /*Used to standardize codon positions. Remember HyPhy indexes from 0, so siteShift = -1 will report first codon as 1*/

/*--Code Overview:--*/
/*
1) Estimate branch lengths by fitting a custom nuc model. These are used throughout.
2) Fit foreground and background codon models, with all parameters except the nonsyn rate tied
	between fg and bg. This is to get estimates of the nucleotide transition rates from a codon model.
3) Site by site meds directional selection analysis: Nuc transition rates are fixed from 2).
	For each site, fit a positive selection null model that allows freely variable syn rate,
	and nonsynBG rates but with nonsynFG constrained to equal syn. Fit a second positive selection
	model, with nonsynFG unconstrained. An LRT is used to determine whether the unconstrained
	positive selection model significantly outperformed the null. So far this is only positive
	selection. The unconstrained positive selection model will be used as the null for tests of
	directional selection. For each of 20 AAs, fit a seperate model that allows a rate multiplier
	favouring substitutions towards that AA. This produces 20 nested models, so a set Bonferroni
	corrected LRTs can be used to identify evidence of direction selection. This batch file outputs
	p-values BEFORE Bonferroni correction.
*/

/*--Code--*/
utility.SetEnvVariable ("ACCEPT_ROOTED_TREES", 1);
utility.SetEnvVariable ("_DO_TREE_REBALANCE_", 0);
utility.SetEnvVariable("COUNT_GAPS_IN_FREQUENCIES", 0);

//ACCEPT_ROOTED_TREES=1;
//_DO_TREE_REBALANCE_ = 0;
//COUNT_GAPS_IN_FREQUENCIES = 0;

//LoadFunctionLibrary ("chooseGeneticCode");


SenseCodons = genetic_code.CountSense(meds.codon_data_info [utility.getGlobalValue("terms.code")]);


/*---------Estimate Branch Lengths Using Nucleotide Model-----------------------*/
//DataSetFilter myFilter = CreateFilter (meds.codon_data,1);
//fprintf(stdout, "%s\n", myFilter);
//HarvestFrequencies (obsNucFreqs, myFilter, 1, 1, 1);





/*---------Begin Setting up custom nuc model---------------*/
/*The parameters that the nucModelString indexes into*/
nucBiasMult = {{"AC*","","AT*","CG*","CT*","GT*"}};

//this pulls out the efv of the gtr results
obsNucFreqs = (meds.gtr_results[utility.getGlobalValue("terms.efv_estimate")])["VALUEINDEXORDER"][0];

/*Initialises and populates a matrix of string values nuc multipliers*/
customRateString = {{"*","","",""}
					{"","*","",""}
					{"","","*",""}
					{"","","","*"}
				 };

for (i=0; i<3; i+=1)
{
	shift = -(i == 0);
	for(j=i+1;j<4;j += 1)
	{
		customRateString[i][j] = Eval ("nucBiasMult["+nucModelString[j+i+shift]+"]");
		customRateString[j][i] = customRateString[i][j];
	}
}



global AC = 1; global AT = 1; global CG = 1; global CT = 1; global GT = 1;

/*To set up a nucleotide model, populate a string version of the rate matrix*/

modelDefString = "";
modelDefString * 16384;

modelDefString* "{" ;
for (i=0; i<4; i += 1)
{
	modelDefString*("{");
	for(j=0;j<4;j=j+1)
	{
		if(j>0)
		{
			modelDefString*(",");
		}
		if(i==j)
		{
			modelDefString*("*");
		}
		else
		{
			modelDefString*(customRateString[i][j]+"t");
		}
	}
	modelDefString*"}";
}
modelDefString*"}";
modelDefString*0;

ExecuteCommands("nucModel = " + modelDefString);
Model RT = (nucModel, obsNucFreqs);

Model FG = (nucModel, obsNucFreqs);
/*This is in case the treestring came with model assignments!*/
/*



//q_ij and Q have to be functions but that seems overly complicated for this????
lfunction FG (type) {
        def = models.DNA.GTR.ModelDescription(type);
        def [utility.getGlobalValue("terms.model.q_ij")] = modelDefString;

	return def;


}
//right now only have a codon data filter and need a nuc one for this but idk how to get the nucleotide one
//model.generic.DefineModel (model_spec, id, arguments, data_filter, estimator_type)
FG = model.generic.DefineModel("models.DNA.GTR.ModelDescription",
        name_space, {
            "0": terms.local
        },
        meds.filter_names,
        None);
*/
//this doesn't work
//model.define_from_components ("FG",modelDefString,obsNucFreqs,1);




/*-----------End of defining nuc model------------*/


/*FIDDLING WITH THE ROOTING. This (commented out) code reroots on an internal
branch and is only certain to work on trees with external branches marked FG.
This is removed and now only rooted trees are accepted.*/

/*
Tree tempTree = treeString;
treeString = RerootTree(treeString,BranchName(tempTree,0));
*/
/*
treeString = trees.GetTreeString(look_for_newick_tree);
Tree givenTree = treeString[^"terms.data.tree"];
*/
fprintf (stdout, "\n\n[PHASE 1. Estimating Branch Lengths using a Nucleotide Model]\n");

/*
LikelihoodFunction theLikFun = (meds.nuc_filter, givenTree, obsFreqs);
Optimize (paramValues, theLikFun);
fprintf (stdout, theLikFun);
*/

meds.partitions_and_trees = trees.LoadAnnotatedTreeTopology.match_partitions(
    meds.codon_data_info[terms.data.partitions],
    meds.name_mapping);

meds.trees = utility.Map(meds.partitions_and_trees, "_partition_", '_partition_[terms.data.tree]');
meds.filter_specification = alignments.DefineFiltersForPartitions(meds.partitions_and_trees, "meds.codon_data" , "meds.filter.", meds.codon_data_info);
meds.filter_names = utility.Map (meds.filter_specification, "_partition_", '_partition_[terms.data.name]');
//meds.store_tree_information();


nucleotide_results = estimators.FitSingleModel_Ext ( meds.filter_names, meds.trees, "FG", meds.gtr_results,
   {
        utility.getGlobalValue  ("terms.run_options.model_type"): utility.getGlobalValue("terms.global"),
        utility.getGlobalValue  ("terms.run_options.retain_lf_object"): TRUE,
        utility.getGlobalValue  ("terms.run_options.retain_model_object"): TRUE
      }
      );


//estimators.BuildLFObject (lf_id, data_filter, tree, model_map, initial_values, model_objects, run_options)
//estimators.BuildLFObject ("busted_sim.lf", busted_sim.filter_names,
//busted_sim.trees, busted_sim.model_map, busted_sim.gtr_results, busted_sim.model_object_map, None);
//meds.model_object_map = { "busted_sim.background" : busted_sim.background.bsrel_model,
                                //    "busted_sim.test" :       busted_sim.test.bsrel_model };

//this doeadn't work yet either
estimators.BuildLFObject("nuc.lf", meds.filter_names, meds.trees,
                          meds.model_map, meds.gtr_results, meds.model_object_map, None);



/*---------end nuc model fit----------*/

/*---------Labels all branches containing branchID with the FG model. This code segment does nothing when the branches are explicitly tagged, and branchMatch doesn't match anything-------------*/
/*--First populate a list of branches containing branchID. branchID could be regExp text--*/
branchMatch = "DONTNAMEANYTAXATHISUNLESSYOUREALLYMEANIT"; /*This is a different (old) way to tag foreground branches - All terminal branches with any substring = branchMatch will get tagged. This is untested and explicit tree tagging is preffered*/
branchID = branchMatch;
BOInames = {};
tips = TipCount(givenTree);
loc = 0;
for (i=0; i<tips; i=i+1)
{
	tempTipName = TipName(givenTree,i);
	tempResult = tempTipName$branchID;
	if (tempResult[1][1]>=0)
		{
			BOInames[loc] = TipName(givenTree,i);
			loc=loc+1;
		}
}
/*--use regular expressions to append the FG model text--*/
numBOI = loc;
newTreeString = treeString;

fprintf (stdout, "The following branches were labeled as foreground: \n");

for (i=0; i< numBOI; i=i+1)
{
	newTreeString = newTreeString^{{BOInames[i]}{BOInames[i]+"{FG}"}};
	/*Just to check what I label*/
	//fprintf (stdout, BOInames[i], "\n");
}
treeString = newTreeString;
/*-----End branch labeling code-----*/

/*---------------------Harvesting Frequencies-----------------------*/
/*---A function that converts a nucFreqMatrix to a vector of freqs--*/
function BuildCodonFrequencies (nucFreqMatrix)
{

	PIStop = 1.0; 		/* denominator */
	result = {SenseCodons,1};    /* resulting codon frequencies */
	hshift = 0;         /* how many stop codons have been counted so far */

	for (h=0; h<64; h=h+1) /* loop over all possible codons */
	{
		first  = h$16;    /* Decompose a codon into 3 nucleotides.
							 The index of the first nucleotide (A=0,C=1,G=2,T=3) is found here,
							 by doing integer division by 16  */
		second = h%16$4;  /* The index of the second nucleotide.
							 First take the remainder of division by 16, i.e. positions 2 and 3
							 and then extract position 2 by integer division by 4*/
		third  = h%4;     /* The index of the third nucleotide.
							 Remainder of integer division by 4*/

						  /* in the end: h = 16*first + 4*second + third */

		if (_Genetic_Code[h]==10) /* stop codon */
		{
			hshift = hshift+1;
			PIStop = PIStop-nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2]; /* adjust the denominator */
		}
		else
		{
			result[h-hshift] = nucFreqMatrix[first][0]*nucFreqMatrix[second][1]*nucFreqMatrix[third][2];
											/* store the frequency for codon h. Notice the substraction of hshift to compensate
											  for the absense of stop codons. The first codon affected by it is
											  TAC (h=49), which gets stored in result[48], because TAA (a stop codon) was skipped. */
		}
	}
	return result*(1.0/PIStop);
}
DataSetFilter 		codonFilter  = CreateFilter	(myData,3,"","","TAA,TAG,TGA"); /* define the codon filter, excluding stop codons */
HarvestFrequencies  (nuc3by4,myData,3,1,1); 			     /* collect position specific nucleotide frequencies */
estimatedCodonFreqs = BuildCodonFrequencies(nuc3by4);
/*----------------------Done with freqs-----------------------------------------*/


/*----------------------Defines a function for creating a custom codon model--------------*/
/*-------Usage: "PopulateModelMatrix ("ModelVarName",Freq3x4,targetAA,customRateString,"nonSynRateTag");"---*/
/*specify targetAA = 21 if you don't want directional selection. NonSynRateTag allows one to control the
name of the "nonsyn" rate variable by appending nonSynRateTag to the end. This allows having seperate foreground
and background nonsyn rates. "customRateString" is NOT a PAUP specifier, but rather a matrix of string valued
multiplers, derived above in the setup of the nuc model*/
/*Part copypasta from MG94customCF3x4.mdl*/

function PopulateModelMatrix (ModelMatrixName&, EFV, targetAA,customRateString,nonSynRateTag)
{
	ModelMatrixDimension = SenseCodons;
	_localNucBiasMult = customRateString;
	ModelMatrixName = {ModelMatrixDimension,ModelMatrixDimension};
	modelDefString = "";
	modelDefString*16384;
	hshift = 0;
	for (h=0; h<64; h=h+1)
	{
		if (_Genetic_Code[h]==10)
		{
			hshift = hshift+1;
			continue;
		}
		vshift = 0;
		for (v = 0; v<64; v=v+1)
		{
			if(h==v)
			{
				continue;
			}
			diff = v-h;
			if (_Genetic_Code[v]==10)
			{
				vshift = vshift+1;
				continue;
			}
			nucPosInCodon = 2;
			if ((h$4==v$4)||((diff%4==0)&&(h$16==v$16))||(diff%16==0)) /* differ by one subsitution only */
			{
				if (h$4==v$4) /* third position */
				{
					transition = v%4;
					transition2= h%4;
				}
				else
				{
					if(diff%16==0) /* first position */
					{
						transition = v$16;
						transition2= h$16;
						nucPosInCodon = 0;
					}
					else   /* second position */
					{
						transition = v%16$4;
						transition2= h%16$4;
						nucPosInCodon = 1;
					}
				}
				hs = Format(h-hshift,0,0);
				vs = Format(v-vshift,0,0);
				ts = Format(transition,0,0);
				ts2= Format(transition2,0,0);
				ps = Format(nucPosInCodon,0,0);
				aa1 = _Genetic_Code[h];
				aa2 = _Genetic_Code[v];

				synOrNon = "nonsyn"+ nonSynRateTag +"*";
				if (aa1==aa2) {synOrNon = "syn*";}
				targetAAmult = "";
				if (aa2==targetAA) {targetAAmult = "((1/(1-omegaT))-1)*";} /*This is a different way of parameterizing the omegaT multiplier that improves optimization speed. Needs to be transformed back at the end*/
				modelDefString*("ModelMatrixName["+hs+"]["+vs+"] := " + targetAAmult + synOrNon + "t*"+_localNucBiasMult[transition][transition2]+"EFV__["+ts+"]["+ps+"];\n");   /*EFV__["+ts+"]["+ps+"] multiplies the 3x4 equilibrium frequency of the column codon*/
			}
	    }
    }
	modelDefString*0;
	ExecuteCommands (modelDefString);
	return 0;
}
/*--------------------End of custom codon model setup function-------------*/

/* --------------------Set up the background model--------------------*/
global syn = 1;
global nonsyn=1;
/*--Populate the transition matrix for the background codon model--*/
PopulateModelMatrix ("MG94xCustomRateMatrix",nuc3by4,21,customRateString,"");
Model MG94xCustom = (MG94xCustomRateMatrix,estimatedCodonFreqs,0);

/* -------------------Set up the foreground model-------------------------*/
global nonsynFG=1;
/*Populate the transition matrix for the foreground codon model*/
PopulateModelMatrix ("MG94xCustomRateMatrixFG",nuc3by4,21,customRateString,"FG");
Model FG = (MG94xCustomRateMatrixFG,estimatedCodonFreqs,0);

/*------------------Assign models to tree------------------------*/
UseModel(MG94xCustom); /*This sets assigns background model to all unlabeled branches*/

/*Just to check my {FG} tags!*/
fprintf (stdout, treeString);

Tree myTreeFG = treeString; /*treeString has "{FG}" appended to all foreground branches*/

/*Forces all branch lengths to be those estimated by the nuc model*/
ReplicateConstraint ("this1.?.t:=this2.?.t__",myTreeFG,givenTree);

/*----------------Optimise likelihood function--------------------*/

fprintf (stdout, "\n\n[PHASE 2. Estimating Substitution Parameters using the Global Codon Model]\n");

LikelihoodFunction lf = (codonFilter, myTreeFG);
Optimize (myRes, lf);
fprintf (stdout, "\n", lf, "\n");

/*--------Forever constrain nuc rates-------*/
AC := AC__; AT := AT__; CG := CG__; CT := CT__; GT := GT__;

/*--------------Setting up the output file--------------------------*/
fprintf (outputFile,CLEAR_FILE,"Site,NullLL,NullBgNonSyn,NullSyn,DivLL,Div_p,DivFgNonSyn,DivBgNonSyn,DivSyn,0LL,0p,0w,0FgNonSyn,0BgNonSyn,0Syn,1LL,1p,1w,1FgNonSyn,1BgNonSyn,1Syn,2LL,2p,2w,2FgNonSyn,2BgNonSyn,2Syn,3LL,3p,3w,3FgNonSyn,3BgNonSyn,3Syn,4LL,4p,4w,4FgNonSyn,4BgNonSyn,4Syn,5LL,5p,5w,5FgNonSyn,5BgNonSyn,5Syn,6LL,6p,6w,6FgNonSyn,6BgNonSyn,6Syn,7LL,7p,7w,7FgNonSyn,7BgNonSyn,7Syn,8LL,8p,8w,8FgNonSyn,8BgNonSyn,8Syn,9LL,9p,9w,9FgNonSyn,9BgNonSyn,9Syn,10LL,10p,10w,10FgNonSyn,10BgNonSyn,10Syn,11LL,11p,11w,11FgNonSyn,11BgNonSyn,11Syn,12LL,12p,12w,12FgNonSyn,12BgNonSyn,12Syn,13LL,13p,13w,13FgNonSyn,13BgNonSyn,13Syn,14LL,14p,14w,14FgNonSyn,14BgNonSyn,14Syn,15LL,15p,15w,15FgNonSyn,15BgNonSyn,15Syn,16LL,16p,16w,16FgNonSyn,16BgNonSyn,16Syn,17LL,17p,17w,17FgNonSyn,17BgNonSyn,17Syn,18LL,18p,18w,18FgNonSyn,18BgNonSyn,18Syn,19LL,19p,19w,19FgNonSyn,19BgNonSyn,19Syn,20LL,20p,20w,20FgNonSyn,20BgNonSyn,20Syn");

/*----------------------------For loop over sites-----------------------------------*/

fprintf (stdout, "\n\n[PHASE 3. Testing for Directional Selection on Foreground Branches]\n");

for(siteIn=1;siteIn<=codonFilter.sites;siteIn += 1)
{
	fprintf (stdout, "Working on site ", siteIn, "\n");

	/*-------------Allow user to select sites and options-------------*/
	site = siteIn +siteShift;
	siteString = "" + (site*3) + "-" + (site*3+2);
	DataSetFilter siteFilter = CreateFilter (myData,3,siteString,"","TAA,TAG,TGA");

	/*-----Count Site Specific Codon and AA Freqs. This is just to exclude certain sites------*/
	HarvestFrequencies (siteCodFreqs, siteFilter, 3, 3, 0);
	AAfreqs = {21,1};
	for (h=0; h<64; h=h+1) /* loop over all possible codons */
	{
		AAfreqs[ _Genetic_Code[h]] += siteCodFreqs[h];
	}

	/*Count how many AA have frequencies greater than 0*/
	numGrtZero = 0;
	for (h=0; h<21; h=h+1) /* loop over all possible codons */
	{
		if(AAfreqs[h]>0)
		{
			numGrtZero = numGrtZero+1;
		}
	}
	/*-------------Only test a site if there is more than one observed AA-------------*/
	if(numGrtZero>1)
	{
		AAlower = 0; AAupper = 20;

		AAlikes = {};
		AAlikesOmegaTs = {};
		AAlikesSyn = {};
		AAlikesBgNonSyn = {};
		AAlikesFgNonSyn = {};

		/*-------------------------loop over amino acid targets--------------------------------*/
		for (AAcount=AAlower; AAcount<AAupper+1; AAcount=AAcount+1)
		{
			/*------Only test for directional selection towards observed AAs-------*/
			if(AAfreqs[AAcount]>0) /*Stop codons are implicitly excluded*/
			{
				fprintf (stdout, "\tTesting target residue ", _hyphyAAOrdering[AAcount], "\n");
				targetAA = AAcount;

				/* --------------------Construct Directional Model---------------------------------*/
				global nonsynDIR=1;
				global omegaT = 0.5; /*reparameterized. 0.5 = 1*/
				PopulateModelMatrix ("MG94xCustomRateMatrixDIR",nuc3by4,targetAA,customRateString,"DIR");
				Model FG = (MG94xCustomRateMatrixDIR,estimatedCodonFreqs,0);

				UseModel(MG94xCustom); /*This assigns background model to all unlabeled branches*/
				Tree myTreeFG = treeString;

				/*Forces all branch lengths to be those estimated by the nuc model*/
				ReplicateConstraint ("this1.?.t:=this2.?.t__",myTreeFG,givenTree);
				omegaT :<1; /*For the reparameterization*/
				syn = 1; /*this is shared between foreground and background*/
				nonsyn=1; /*background nonsyn rate*/

				LikelihoodFunction lf = (siteFilter, myTreeFG);
				Optimize (mySiteRes, lf);
				/*fprintf (stdout, "\n", lf, "\n");*/

				unConstrainedRatio = mySiteRes[1][0];
				AAlikes[AAcount] = unConstrainedRatio;
				AAlikesOmegaTs[AAcount] = omegaT;
				AAlikesSyn[AAcount] = syn;
				AAlikesBgNonSyn[AAcount] = nonsyn;
				AAlikesFgNonSyn[AAcount] = nonsynDIR;
			}
			else
			{
				fprintf (stdout, "\tSkipping unobserved target residue ", _hyphyAAOrdering[AAcount], "\n");
				AAlikes[AAcount] = -99999999;
				AAlikesOmegaTs[AAcount] = -99;
				AAlikesSyn[AAcount] = -99;
				AAlikesBgNonSyn[AAcount] = -99;
				AAlikesFgNonSyn[AAcount] = -99;
			}
		}

		/*---------------------Set up model allowing non-neutral selection on FG---------------*/
		syn = 1;
		nonsyn=1;
		nonsynDIR=1; /*Allow positive selection on FG*/
		omegaT :=0.5; /*Constrain omegaT to 1 for null model - reparameterized*/
		LikelihoodFunction lf = (siteFilter, myTreeFG);
		Optimize (mySiteRes, lf);
		constrainedRatioPos = mySiteRes[1][0];
		DivFgNonSyn = nonsynDIR;
		DivBgNonSyn = nonsyn;
		DivSyn = syn;
		/*--End non-neutral model--*/

		/*---------------------Set up null model forcing neutral selection on FG----------------*/
		syn = 1;
		nonsyn=1;
		nonsynDIR:=syn; /*Force neutral selection on FG*/
		omegaT :=0.5; /*Constrain omegaT to 1 for null model - reparameterized*/
		LikelihoodFunction lf = (siteFilter, myTreeFG);
		Optimize (mySiteRes, lf);
		/*fprintf (stdout, "\n", lf, "\n");*/
		constrainedRatioNoPos = mySiteRes[1][0];
		/*--End null model--*/

		/*-----------------------Display and Write Results------------------------*/
		/*---"Site,NullLL,NullBgNonSyn,NullSyn,DivLL,Div_p,DivFgNonSyn,DivBgNonSyn,DivSyn,0LL,0p,0w,0FgNonSyn,0BgNonSyn,0Syn"---*/
		outputString = ""+siteIn;
		outputString = outputString + "," + constrainedRatioNoPos + "," + nonsyn + "," + syn;
		fprintf (stdout, "\nTests for site ", siteIn,"\n");

		/*--First a test for general positive selection at the site--*/
		lrtPosVsNoPos = 2*(constrainedRatioPos-constrainedRatioNoPos);
		pValPosVsNoPos = 1-CChi2 (lrtPosVsNoPos, 1);
		fprintf (stdout, "\tLikelihood Ratio Test for diversifying selection");
		fprintf (stdout, ": ", lrtPosVsNoPos);
		fprintf (stdout, "  p-value: ", pValPosVsNoPos, "\n");

		outputString = outputString + "," + constrainedRatioNoPos + "," + pValPosVsNoPos + "," + DivFgNonSyn + "," + DivBgNonSyn + "," + DivSyn;

		/*--Now test for directional selection vs positive selection--*/
		fprintf (stdout, "\tTesting for directional selection againts null allowing positive selection in foreground", "\n");
		for (AAcount=AAlower; AAcount<AAupper+1; AAcount=AAcount+1)
		{
			lrtScore = 2*(AAlikes[AAcount]-constrainedRatioPos);
			pValues = 1-CChi2 (2*(AAlikes[AAcount]-constrainedRatioPos), 1);
			if (AAfreqs[AAcount] > 0)
			{
				fprintf (stdout, "\tLikelihood Ratio Test for: ", _hyphyAAOrdering[AAcount]);
				fprintf (stdout, "  : ", lrtScore);
				fprintf (stdout, "  p-value: ", pValues, "\n");
			}
			outputString = outputString + "," + AAlikes[AAcount] + "," + pValues + "," + AAlikesOmegaTs[AAcount] + "," + AAlikesFgNonSyn[AAcount] + "," + AAlikesBgNonSyn[AAcount] + "," + AAlikesSyn[AAcount];
		}
		/*Testing for directional vs neutral selection can be done in post-processesing*/
		fprintf (outputFile,"\n",outputString);
		/*Clear constraints to test another site*/
		ClearConstraints(omegaT);
		ClearConstraints(nonsynDIR);
	}
	else
	{
		fprintf (stdout,"Skipped site ",siteIn," because it is invariable\n");
	}

} /*End main loop*/
