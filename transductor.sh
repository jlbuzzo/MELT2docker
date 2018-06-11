#!/bin/bash


# Arguments validation.
if (( $# < 4 )); then
  	echo -e "Error in $0: Not enough arguments." >&2
	exit 1;
fi

# Reference paths.
CONFIG_FILE=config.yml

BASE=/home
IMG_CONFIG=$BASE/config
IMG_IN=$BASE/inputs
IMG_OUT=$BASE/outputs
IMG_REF=$BASE/reference
IMG_BED=$BASE/MELTv2.1.4/add_bed_files/Hg38/Hg38.genes.bed

declare -A PATH_H DIR_H FILE_H

while getopts ":c:i:o:r:" OPTION; do
	case "$OPTION" in
		c) PATH_H["c"]=$OPTARG ;;
		i) PATH_H["i"]=$OPTARG ;;
		o) PATH_H["o"]=$OPTARG ;;
		r) PATH_H["r"]=$OPTARG ;;
	esac
done

for var in $(echo -e "c\ti\to\tr"); do
	if [[ ! -v PATH_H[$var] ]]; then
		echo -e "No flag '$var' defined." >&2
		exit 1
	fi
done

for var in $(echo -e "c\ti\to\tr"); do
	path="${PATH_H[$var]}"
	DIR_H[$var]="${path%/*}"
	FILE_H[$var]="${path##*/}"
	
	# Debug code.
	#echo "dir ${DIR_H[$var]} file ${FILE_H[$var]}"
done

cat << EOF > $CONFIG_FILE
EOF

# Run docker.
docker run \
	--rm -ti \
	-u $(id -u):$(id -g) \
	-v "${DIR_H["c"]}:$IMG_CONFIG" \
	-v "${DIR_H["i"]}:$IMG_IN" \
	-v "${DIR_H["o"]}:$IMG_OUT" \
	-v "${DIR_H["r"]}:$IMG_REF" \
	melt
#EOF
