#!/bin/bash


while [ $# -gt 0 ]
  do
           case "$1" in
                  --cov) cov="$2"; shift;;
                  --manifest-summary) manifest_summary="$2"; shift;;
                  --competitive-mapping-json) competitive_mapping_json="$2"; shift;;
                  --) shift;;
           esac
           shift;
  done

# Generate JSON output
csvjoin --columns "#rname" "${cov}" "${manifest_summary}" > joined.csv

# rename endpos to length (this will occur throughout the file, so could be improved)
sed -i '1s/endpos/length/' joined.csv

cat joined.csv |
csvcut --delimiter="," -c "#rname","genome_name","length","coverage","numreads","meandepth" |
csvsort -r -c "coverage" |
csvjson -i 4 > "${competitive_mapping_json}"
