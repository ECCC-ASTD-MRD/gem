### Mandatory checklist:

<!--For each point below, choose 'YES' or 'NO' -->

* new functionality?  **YES/NO**
* fix bug in existing functionality? **YES/NO**
* changes to namelist variables? (addition/removal/modification) **YES/NO**
* changes to input and output files? (new file/filename change/removal of file) **YES/NO**
* change to results (contents of output files)? **YES/NO**
* changes respect the coding standards (follow link below to "Coding standards Top 10")? **YES/NO**

    * https://wiki.cmc.ec.gc.ca/wiki/MIDAS/Coding_Standards_Top_10

### Description of changes:

<!--The text here should describe how the change was implemented.-->
<!--Detail here the changes answered as YES in the previous section-->

### Description of the impact on the results

<!--Describe how the results are affected by the code modifications introduced
in this merge request:

 1. quick explanation of why the results are affected:
   * e.g. order of obs changed affecting cost function calculation
 2. list the program(s) that are affected
   * e.g. `midas-var.Abs`
 3. describe in which configurations the impacts are seen
   * e.g. all operational NWP systems using `midas-var.Abs` program
 4. level of significance of the changes to the results (choose one of following)
   * no impact (completely backward compatible)
   * minor, only due to numerical round-off error or
   * major, impacting the meteorological evaluation
      * give the link to results of data assimilation experiments in
        which the impact of the changes were evaluated
-->

Delete this line and replace it with your description of the impact on
the results based on the text above. If there are no absolutely no impacts
on the results then simply state "No impact"

### Addition to CHANGELOG:

<!--Some oneliners describing changes for the whole merge-request-->
<!--That information will be added to the 'CHANGELOG.md' file-->
<!--Put any information relevant to the user, especially non-backward compatible changes-->
<!--   * new functionality  -->
<!--   * Namelist variables -->
<!--   * input/output files -->
<!--   * results            -->

Delete this line and replace it with the CHANGELOG entry that you put
in the [CHANGELOG under one of the subheadings after "Unreleased"](CHANGELOG.md#Unreleased).

If the changes have any impact on the results, then you must include either
  * "minor impact on results" or
  * "major impact on results".

### Resolved issues:

<!--Put the list of issues that this merge request resolves-->
Closes #...

<!--(choose one of the following labels)-->
/label ~Optimization
/label ~Bug
/label ~FeatureRequest