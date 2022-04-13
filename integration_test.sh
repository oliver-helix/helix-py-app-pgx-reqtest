#!/bin/bash

set -euo pipefail

BIN_DIR=$(dirname $0)

# Parse command line options.
function usage {
    local exit_code=$1
    echo "Usage: $0 [-h] [-c] [arg [arg ...]]

            -h  Show this help message and exit.
            -c  Collect code coverage data and write results to htmlcov/.
           arg  Extra positional argument passed to call_pgx to override config settings.
             " >&2
    exit $exit_code
}
GET_COVERAGE=''  # By default, do not collect code coverage.
while getopts :hc option; do
    case "${option}" in
        h)
            usage 0
            ;;
        c)
            GET_COVERAGE=1
            ;;
        :)
            echo "Argument required for option -$OPTARG" >&2
            usage 1
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage 1
            ;;
    esac
done
shift $((OPTIND - 1))  # Remove parsed args from $@, leaving only any positional args.

set -x
for ASSAY_VERSION in v2 v3 v4; do

    TEST_COV_DIR=htmlcov_$ASSAY_VERSION
    METADATA='s3://helix-biofx-data-batch/helix-py-app-pgx/integration_test/integration_test_'$ASSAY_VERSION'.json'
    observed_dir='/integration_test/'$ASSAY_VERSION'/output'

    echo -e "\n== Running integration test for $ASSAY_VERSION =="
    cd $BIN_DIR
    if [ $GET_COVERAGE ]; then
        CALL_PGX="coverage run $(which call_pgx)"
        rm -rf $TEST_COV_DIR
    else
        CALL_PGX='call_pgx'
    fi

    $CALL_PGX \
        --input_json $METADATA \
        --platform 'batch' \
        --output_json 'integration_test_directory.json' \
        --s3_dir 'test_s3' \
        --genes 'CYP2D6,CYP2C19,VKORC1' \
        app.skip_uploading=yes \
        parameters.low_baseline_size=10 \
        processes.download=4 \
        processes.counts_data=4 \
        processes.variant_calling=4 \
        output.directory=$observed_dir \
        "$@"
    if [ $GET_COVERAGE ]; then
        coverage combine
        coverage html --directory=$TEST_COV_DIR
    fi
    set +x

    echo -e "\n== Checking output files against expected files... =="
    expected_dir=$BIN_DIR/integration_test_expected/$ASSAY_VERSION
    pass=1
    while read expected_path; do
        name=${expected_path#$expected_dir/}
        observed_path=$observed_dir/$name

        file_pass=1
        if [[ $name =~ \.json$ ]]; then
            # Tolerate different pgxVersion for *.json files.
            diff -q <(grep -v '"pgxVersion"' $expected_path) \
                    <(grep -v '"pgxVersion"' $observed_path) \
                || { file_pass=0; echo "Files $expected_path and $observed_path differ materially"; }
        else
            # Other files must match exactly.
            diff -q $expected_path $observed_path || file_pass=0
        fi

        if [ $file_pass -eq 1 ]; then
            printf '%-65s  OK\n' $name
        else
            # Failure for any file results in failure of the whole integration test.
            pass=0
        fi
    done < <(find $expected_dir -type f | sort)

    if [ $pass -eq 1 ]; then
        echo -e "\n== Integration test for $ASSAY_VERSION passed! =="
    else
        echo -e "\n== Integration test $ASSAY_VERSION failed. =="
        exit 1
    fi
done
