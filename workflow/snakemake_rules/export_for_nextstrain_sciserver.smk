# How to run: if no region is specified, it'll run a subsampled global build (120 per division)
# If a region is selected, it'll do 280/division for that region, and 20/division in the rest of the world
#       -- preferentially sequences near the focal sequences
#
# To run a regional build, be sure to update the list of regions in `config/nextstrain_profiles.yaml`.
#
# You can run all builds in parallel!
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
#
# Or you can specify final or intermediate output files like so:
#   snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_europe.json (subsampled regional focal)
#   snakemake --profile nextstrain_profiles/nextstrain auspice/ncov_global.json (subsampled global)
#
# To update ordering/lat_longs after AWS download:
#   snakemake --touch --forceall --profile nextstrain_profiles/nextstrain
#   snakemake --profile nextstrain_profiles/nextstrain clean_export_regions
#   snakemake --profile nextstrain_profiles/nextstrain export_all_regions
# When done adjusting lat-longs & orders, remember to run
#   snakemake --profile nextstrain_profiles/nextstrain all_regions
# to produce the final Auspice files!

def get_todays_date():
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d')
    return date


rule all_regions:
    input:
        auspice_json = expand(OUTDIR + "/auspice/ncov_{build_name}.json", build_name=BUILD_NAMES),
        tip_frequencies_json = expand(OUTDIR + "/auspice/ncov_{build_name}_tip-frequencies.json", build_name=BUILD_NAMES),
        dated_auspice_json = expand(OUTDIR + "/auspice/ncov_{build_name}_{date}.json", build_name=BUILD_NAMES, date=get_todays_date()),
        dated_tip_frequencies_json = expand(OUTDIR + "/auspice/ncov_{build_name}_{date}_tip-frequencies.json", build_name=BUILD_NAMES, date=get_todays_date()),
        staged_dated_auspice_json = expand(STAGE_DIR + "/ncov-" + STAGE_NAME + "_" + RUN_NAME + "_{build_name}_{date}.json", build_name=BUILD_NAMES, date=get_todays_date()),
        staged_dated_tip_frequencies_json = expand(STAGE_DIR + "/ncov-" + STAGE_NAME + "_" + RUN_NAME + "_{build_name}_{date}_tip-frequencies.json", build_name=BUILD_NAMES, date=get_todays_date())        
        
        #staged_dated_auspice_json = expand(STAGE_DIR + "/ncov_{build_name}_{date}.json", build_name=BUILD_NAMES, date=get_todays_date()),
        #staged_dated_tip_frequencies_json = expand(STAGE_DIR + "/ncov_{build_name}_{date}_tip-frequencies.json", build_name=BUILD_NAMES, date=get_todays_date())        

# This cleans out files to allow re-run of 'normal' run with `export`
# to check lat-longs & orderings
rule clean_export_regions:
    message: "Removing export files: {input}"
    params:
        *expand(OUTDIR + "/results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        *expand(OUTDIR + "/results/{build_name}/colors.tsv", build_name=BUILD_NAMES)
    conda: config["conda_environment"]
    shell:
        "rm -f {params}"

# Allows 'normal' run of export to be forced to correct lat-long & ordering
# Runs an additional script to give a list of locations that need colors and/or lat-longs
rule export_all_regions:
    input:
        auspice_json = expand(OUTDIR + "/results/{build_name}/ncov_with_accessions.json", build_name=BUILD_NAMES),
        lat_longs = config["files"]["lat_longs"],
        metadata = [_get_metadata_by_build_name(build_name).format(build_name=build_name)
                    for build_name in BUILD_NAMES],
        colors = expand(OUTDIR + "/results/{build_name}/colors.tsv", build_name=BUILD_NAMES),
    conda: config["conda_environment"]
    shell:
        """
        python3 ./scripts/check_missing_locations.py \
            --metadata {input.metadata} \
            --colors {input.colors} \
            --latlong {input.lat_longs}
        """

rule all_mutation_frequencies:
    input: expand(OUTDIR + "/results/{build_name}/nucleotide_mutation_frequencies.json", build_name=BUILD_NAMES)

#
# Rules for custom auspice exports for the Nextstrain team.
#

rule dated_json:
    message: "Copying dated Auspice JSON"
    input:
        auspice_json = rules.finalize.output.auspice_json,
        tip_frequencies_json = rules.tip_frequencies.output.tip_frequencies_json
    output:
        dated_auspice_json = OUTDIR + "/auspice/ncov_{build_name}_{date}.json",
        dated_tip_frequencies_json = OUTDIR + "/auspice/ncov_{build_name}_{date}_tip-frequencies.json",
        staged_dated_auspice_json = STAGE_DIR + "/ncov-" + STAGE_NAME + "_" + RUN_NAME + "_{build_name}_{date}.json",
        staged_dated_tip_frequencies_json = STAGE_DIR + "/ncov-" + STAGE_NAME + "_" + RUN_NAME + "_{build_name}_{date}_tip-frequencies.json"       
    conda: config["conda_environment"]
    shell:
        """
        cp {input.auspice_json} {output.dated_auspice_json}
        cp {input.tip_frequencies_json} {output.dated_tip_frequencies_json}
        cp {input.auspice_json} {output.staged_dated_auspice_json}
        cp {input.tip_frequencies_json} {output.staged_dated_tip_frequencies_json}
        """

#
# Deployment and error handlers, including Slack messaging integrations.
#

from os import environ

#curl -X POST http://172.23.38.10:9001/restart

rule deploy_to_staging:
    input:
        *rules.all_regions.input
    params:
        staging_message = "Deployed to https://alpha.nextstrain.idies.jhu.edu",
        staging_url = config["staging_url"]
    output:
        deploy_report = OUTDIR + "/results/deploy_report.txt"
    conda: config["conda_environment"]
    shell:
        """
         curl -X POST "{params.staging_url}/restart"  
         echo {params.staging_message} > {output.deploy_report}
        """
